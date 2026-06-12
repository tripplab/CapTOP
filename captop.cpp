// captop.cpp
// Build:
//   g++ -std=c++17 -O2 -Wall -Wextra -pedantic captop.cpp -o captop
//
// Usage:
//   ./captop --version
//   ./captop validate mesh.msh
//   ./captop validate mesh.msh --grid strict-cube
//   ./captop validate mesh.msh --grid approximate-cube
//   ./captop validate mesh.msh --grid rectilinear
//   ./captop validate mesh.msh --tol 1e-8
//   ./captop validate mesh.msh --cube-rel-tol 0.02
//   ./captop validate mesh.msh --ignore-non-hexa

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace captop {

static const char* CAPTOP_VERSION = "0.1.1-stage3";

enum class GridMode {
    StrictCube,
    ApproximateCube,
    Rectilinear
};

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Node {
    long long original_id = 0;
    Vec3 p;
};

struct HexElement {
    long long original_id = 0;
    std::array<long long, 8> node_ids{};
    bool has_material = false;
    long long material = 0;
    std::size_t mesh_block_id = 0;
    std::size_t line_number = 0;
};

struct MeshBlock {
    std::size_t id = 0;
    std::size_t line_number = 0;
    int dimension = 3;
    std::string elemtype;
    int nnode = 0;
    bool supported_hexa = false;
    std::size_t element_count = 0;
};

struct Mesh {
    std::map<long long, Node> nodes;
    std::vector<HexElement> hexes;
    std::vector<MeshBlock> blocks;
    std::vector<std::string> warnings;
};

struct ParseOptions {
    bool ignore_non_hexa = false;
};

struct ValidateOptions {
    GridMode grid_mode = GridMode::StrictCube;
    double tol = 1.0e-8;
    double cube_rel_tol = 0.02;
    std::size_t max_errors = 30;
};

struct ScalarStats {
    std::size_t count = 0;
    double min = 0.0;
    double max = 0.0;
    double mean = 0.0;
    double stddev = 0.0;
    double q1 = 0.0;
    double median = 0.0;
    double q3 = 0.0;
};

struct WorstCell {
    long long element_id = 0;
    std::size_t line_number = 0;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double s_min = 0.0;
    double s_max = 0.0;
    double rel_spread = 0.0;
};

struct SideLengthStats {
    bool available = false;
    ScalarStats dx;
    ScalarStats dy;
    ScalarStats dz;
    ScalarStats rel_spread;
    std::size_t above_0_1_percent = 0;
    std::size_t above_0_5_percent = 0;
    std::size_t above_1_0_percent = 0;
    std::size_t above_2_0_percent = 0;
    std::size_t above_5_0_percent = 0;
    double max_abs_spread = 0.0;
    double max_abs_dx_dy = 0.0;
    double max_abs_dx_dz = 0.0;
    double max_abs_dy_dz = 0.0;
    double mean_abs_dx_dy = 0.0;
    double mean_abs_dx_dz = 0.0;
    double mean_abs_dy_dz = 0.0;
    std::vector<WorstCell> worst_cells;
};

struct AxisIntervalStats {
    bool available = false;
    std::size_t coordinates = 0;
    std::size_t intervals = 0;
    double min_interval = 0.0;
    double max_interval = 0.0;
    double mean_interval = 0.0;
    double median_interval = 0.0;
    std::size_t unique_intervals = 0;
};

struct RectilinearAxisStats {
    bool available = false;
    AxisIntervalStats x;
    AxisIntervalStats y;
    AxisIntervalStats z;
    Vec3 median_interval;
};

struct ValidationResult {
    bool valid = false;

    std::vector<std::string> errors;
    std::vector<std::string> warnings;

    std::size_t n_nodes = 0;
    std::size_t n_hexes = 0;
    std::size_t n_blocks = 0;
    std::size_t n_supported_blocks = 0;
    std::size_t n_unsupported_blocks = 0;

    std::size_t occupied_cubes = 0;
    std::size_t face_adjacencies = 0;
    std::size_t boundary_faces = 0;

    long long nx = 0;
    long long ny = 0;
    long long nz = 0;

    Vec3 min_bounds;
    Vec3 max_bounds;
    Vec3 origin;
    Vec3 spacing;

    bool has_bounds = false;

    SideLengthStats side_stats;
    RectilinearAxisStats axis_stats;
    bool topology_compatible = false;

    std::map<long long, std::size_t> material_counts;
};

struct ParseError : public std::runtime_error {
    explicit ParseError(const std::string& msg) : std::runtime_error(msg) {}
};

struct ValidationError : public std::runtime_error {
    explicit ValidationError(const std::string& msg) : std::runtime_error(msg) {}
};

static std::string trim(const std::string& s) {
    const std::string ws = " \t\n\r";
    const std::size_t a = s.find_first_not_of(ws);
    if (a == std::string::npos) {
        return "";
    }
    const std::size_t b = s.find_last_not_of(ws);
    return s.substr(a, b - a + 1);
}

static std::string to_lower(std::string s) {
    std::transform(
        s.begin(),
        s.end(),
        s.begin(),
        [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        }
    );
    return s;
}

static bool starts_with_case_insensitive(const std::string& s, const std::string& prefix) {
    if (s.size() < prefix.size()) {
        return false;
    }
    return to_lower(s.substr(0, prefix.size())) == to_lower(prefix);
}

static std::vector<std::string> tokenize_gid_line(const std::string& line) {
    std::string expanded;
    expanded.reserve(line.size() + 8);

    for (char c : line) {
        if (c == '=' || c == ',' || c == ';') {
            expanded.push_back(' ');
            expanded.push_back(c);
            expanded.push_back(' ');
        } else {
            expanded.push_back(c);
        }
    }

    std::istringstream iss(expanded);
    std::vector<std::string> tokens;
    std::string tok;

    while (iss >> tok) {
        if (tok == "=" || tok == "," || tok == ";") {
            continue;
        }
        tokens.push_back(tok);
    }

    return tokens;
}

static bool parse_long_long(const std::string& s, long long& out) {
    char* end = nullptr;
    errno = 0;
    const long long v = std::strtoll(s.c_str(), &end, 10);
    if (errno != 0 || end == s.c_str() || *end != '\0') {
        return false;
    }
    out = v;
    return true;
}

static bool parse_int(const std::string& s, int& out) {
    long long tmp = 0;
    if (!parse_long_long(s, tmp)) {
        return false;
    }
    if (tmp < static_cast<long long>(std::numeric_limits<int>::min()) ||
        tmp > static_cast<long long>(std::numeric_limits<int>::max())) {
        return false;
    }
    out = static_cast<int>(tmp);
    return true;
}

static bool parse_double(const std::string& s, double& out) {
    char* end = nullptr;
    errno = 0;
    const double v = std::strtod(s.c_str(), &end);
    if (errno != 0 || end == s.c_str() || *end != '\0' || !std::isfinite(v)) {
        return false;
    }
    out = v;
    return true;
}

static bool nearly_equal(double a, double b, double tol) {
    const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
    return std::abs(a - b) <= tol * scale;
}

static bool positive_length(double x, double tol) {
    return x > tol * std::max(1.0, std::abs(x));
}

static std::string line_context(std::size_t line_number, const std::string& msg) {
    std::ostringstream oss;
    oss << "line " << line_number << ": " << msg;
    return oss.str();
}

static MeshBlock parse_mesh_header(
    const std::string& line,
    std::size_t line_number,
    std::size_t block_id
) {
    const std::vector<std::string> tokens = tokenize_gid_line(line);

    MeshBlock block;
    block.id = block_id;
    block.line_number = line_number;
    block.dimension = 3;
    block.elemtype = "";
    block.nnode = 0;

    for (std::size_t i = 0; i < tokens.size(); ++i) {
        const std::string key = to_lower(tokens[i]);

        if (key == "dimension") {
            if (i + 1 >= tokens.size()) {
                throw ParseError(line_context(line_number, "missing value after 'dimension'"));
            }
            int dim = 0;
            if (!parse_int(tokens[i + 1], dim)) {
                throw ParseError(line_context(line_number, "invalid dimension value '" + tokens[i + 1] + "'"));
            }
            block.dimension = dim;
            ++i;
        } else if (key == "elemtype") {
            if (i + 1 >= tokens.size()) {
                throw ParseError(line_context(line_number, "missing value after 'elemtype'"));
            }
            block.elemtype = tokens[i + 1];
            ++i;
        } else if (key == "nnode") {
            if (i + 1 >= tokens.size()) {
                throw ParseError(line_context(line_number, "missing value after 'nnode'"));
            }
            int nn = 0;
            if (!parse_int(tokens[i + 1], nn)) {
                throw ParseError(line_context(line_number, "invalid nnode value '" + tokens[i + 1] + "'"));
            }
            block.nnode = nn;
            ++i;
        }
    }

    const std::string et = to_lower(block.elemtype);
    block.supported_hexa = (block.dimension == 3 && et == "hexahedra" && block.nnode == 8);

    return block;
}

static Mesh parse_gid_mesh(const std::string& path, const ParseOptions& opts) {
    std::ifstream in(path.c_str());
    if (!in) {
        throw ParseError("cannot open input file '" + path + "'");
    }

    enum class State {
        Outside,
        Coordinates,
        Elements
    };

    Mesh mesh;
    State state = State::Outside;
    std::size_t current_block_index = std::numeric_limits<std::size_t>::max();

    std::string line;
    std::size_t line_number = 0;

    while (std::getline(in, line)) {
        ++line_number;

        const std::string stripped = trim(line);
        if (stripped.empty()) {
            continue;
        }

        if (starts_with_case_insensitive(stripped, "#")) {
            continue;
        }

        const std::string lower = to_lower(stripped);

        if (starts_with_case_insensitive(lower, "mesh")) {
            MeshBlock block = parse_mesh_header(stripped, line_number, mesh.blocks.size());
            mesh.blocks.push_back(block);
            current_block_index = mesh.blocks.size() - 1;
            state = State::Outside;

            if (!block.supported_hexa) {
                std::ostringstream oss;
                oss << "line " << line_number
                    << ": mesh block " << block.id
                    << " is not supported as a cubic/hexahedral block"
                    << " (dimension=" << block.dimension
                    << ", elemtype='" << block.elemtype
                    << "', nnode=" << block.nnode << ")";
                mesh.warnings.push_back(oss.str());
            }

            continue;
        }

        if (lower == "coordinates" || starts_with_case_insensitive(lower, "coordinates ")) {
            if (current_block_index == std::numeric_limits<std::size_t>::max()) {
                throw ParseError(line_context(line_number, "'coordinates' found before any 'mesh' header"));
            }
            state = State::Coordinates;
            continue;
        }

        if (lower == "end coordinates" || lower == "endcoordinates") {
            state = State::Outside;
            continue;
        }

        if (lower == "elements" || starts_with_case_insensitive(lower, "elements ")) {
            if (current_block_index == std::numeric_limits<std::size_t>::max()) {
                throw ParseError(line_context(line_number, "'elements' found before any 'mesh' header"));
            }
            state = State::Elements;
            continue;
        }

        if (lower == "end elements" || lower == "endelements") {
            state = State::Outside;
            continue;
        }

        if (state == State::Coordinates) {
            const MeshBlock& block = mesh.blocks.at(current_block_index);
            const std::vector<std::string> tokens = tokenize_gid_line(stripped);

            if (tokens.size() < 1 + static_cast<std::size_t>(block.dimension)) {
                throw ParseError(line_context(line_number, "malformed coordinate record"));
            }

            long long node_id = 0;
            if (!parse_long_long(tokens[0], node_id)) {
                throw ParseError(line_context(line_number, "invalid node id '" + tokens[0] + "'"));
            }

            if (block.dimension != 3) {
                throw ParseError(line_context(line_number, "only dimension 3 coordinates are supported"));
            }

            double x = 0.0;
            double y = 0.0;
            double z = 0.0;

            if (!parse_double(tokens[1], x) ||
                !parse_double(tokens[2], y) ||
                !parse_double(tokens[3], z)) {
                throw ParseError(line_context(line_number, "invalid coordinate value"));
            }

            Node n;
            n.original_id = node_id;
            n.p = Vec3{x, y, z};

            const auto it = mesh.nodes.find(node_id);
            if (it != mesh.nodes.end()) {
                const Vec3& old = it->second.p;
                if (old.x != x || old.y != y || old.z != z) {
                    throw ParseError(
                        line_context(
                            line_number,
                            "node id " + std::to_string(node_id) + " is defined more than once with different coordinates"
                        )
                    );
                }
                mesh.warnings.push_back(
                    line_context(
                        line_number,
                        "node id " + std::to_string(node_id) + " is repeated with identical coordinates; keeping first definition"
                    )
                );
            } else {
                mesh.nodes.insert(std::make_pair(node_id, n));
            }

            continue;
        }

        if (state == State::Elements) {
            MeshBlock& block = mesh.blocks.at(current_block_index);
            const std::vector<std::string> tokens = tokenize_gid_line(stripped);

            if (tokens.empty()) {
                continue;
            }

            block.element_count += 1;

            if (!block.supported_hexa) {
                if (opts.ignore_non_hexa) {
                    continue;
                }

                std::ostringstream oss;
                oss << "element record belongs to unsupported block " << block.id
                    << " (dimension=" << block.dimension
                    << ", elemtype='" << block.elemtype
                    << "', nnode=" << block.nnode << ")";
                throw ParseError(line_context(line_number, oss.str()));
            }

            if (tokens.size() < 9) {
                throw ParseError(line_context(line_number, "hexahedral element requires id plus 8 node ids"));
            }

            HexElement h;
            h.mesh_block_id = current_block_index;
            h.line_number = line_number;

            if (!parse_long_long(tokens[0], h.original_id)) {
                throw ParseError(line_context(line_number, "invalid element id '" + tokens[0] + "'"));
            }

            for (std::size_t k = 0; k < 8; ++k) {
                long long nid = 0;
                if (!parse_long_long(tokens[1 + k], nid)) {
                    throw ParseError(line_context(line_number, "invalid node id '" + tokens[1 + k] + "' in element"));
                }
                h.node_ids[k] = nid;
            }

            if (tokens.size() >= 10) {
                long long mat = 0;
                if (!parse_long_long(tokens[9], mat)) {
                    throw ParseError(line_context(line_number, "invalid material/layer id '" + tokens[9] + "'"));
                }
                h.has_material = true;
                h.material = mat;
            }

            if (tokens.size() > 10) {
                mesh.warnings.push_back(
                    line_context(
                        line_number,
                        "extra columns after material/layer id are ignored"
                    )
                );
            }

            mesh.hexes.push_back(h);
            continue;
        }

        mesh.warnings.push_back(line_context(line_number, "ignored line outside recognized GiD sections"));
    }

    return mesh;
}

static std::vector<double> unique_sorted_tol(std::vector<double> values, double tol) {
    std::sort(values.begin(), values.end());

    std::vector<double> unique;
    for (double v : values) {
        if (unique.empty()) {
            unique.push_back(v);
        } else if (!nearly_equal(unique.back(), v, tol)) {
            unique.push_back(v);
        }
    }

    return unique;
}

static int classify_to_pair(double v, double a, double b, double tol) {
    const bool near_a = nearly_equal(v, a, tol);
    const bool near_b = nearly_equal(v, b, tol);

    if (near_a && !near_b) {
        return 0;
    }

    if (near_b && !near_a) {
        return 1;
    }

    if (near_a && near_b) {
        return 0;
    }

    return -1;
}

struct CellGeom {
    long long element_id = 0;
    std::size_t line_number = 0;
    bool has_material = false;
    long long material = 0;

    double xmin = 0.0;
    double xmax = 0.0;
    double ymin = 0.0;
    double ymax = 0.0;
    double zmin = 0.0;
    double zmax = 0.0;

    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
};


struct CellIndex {
    long long i = 0;
    long long j = 0;
    long long k = 0;
};

struct CellIndexHash {
    std::size_t operator()(const CellIndex& c) const {
        const std::size_t h1 = std::hash<long long>()(c.i);
        const std::size_t h2 = std::hash<long long>()(c.j);
        const std::size_t h3 = std::hash<long long>()(c.k);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6U) + (h1 >> 2U))
                  ^ (h3 + 0x9e3779b97f4a7c15ULL + (h2 << 6U) + (h2 >> 2U));
    }
};

struct CellIndexEqual {
    bool operator()(const CellIndex& a, const CellIndex& b) const {
        return a.i == b.i && a.j == b.j && a.k == b.k;
    }
};

static std::string cell_key_string(const CellIndex& c) {
    std::ostringstream oss;
    oss << "(" << c.i << "," << c.j << "," << c.k << ")";
    return oss.str();
}

static void add_validation_error(
    ValidationResult& result,
    const std::string& msg,
    std::size_t max_errors
) {
    if (result.errors.size() < max_errors) {
        result.errors.push_back(msg);
    } else if (result.errors.size() == max_errors) {
        result.errors.push_back("maximum error count reached; additional errors suppressed");
    }
}

static bool coordinate_to_strict_index(
    double value,
    double origin,
    double h,
    double tol,
    long long& out
) {
    const double raw = (value - origin) / h;
    const double rounded = std::round(raw);

    if (rounded < static_cast<double>(std::numeric_limits<long long>::min()) ||
        rounded > static_cast<double>(std::numeric_limits<long long>::max())) {
        return false;
    }

    const double expected = origin + rounded * h;
    if (!nearly_equal(value, expected, tol)) {
        return false;
    }

    out = static_cast<long long>(rounded);
    return true;
}

static bool find_axis_index(
    const std::vector<double>& axis,
    double value,
    double tol,
    long long& out
) {
    if (axis.empty()) {
        return false;
    }

    auto it = std::lower_bound(axis.begin(), axis.end(), value);

    std::size_t best = 0;
    double best_diff = std::numeric_limits<double>::infinity();

    if (it != axis.end()) {
        const std::size_t idx = static_cast<std::size_t>(std::distance(axis.begin(), it));
        best = idx;
        best_diff = std::abs(axis[idx] - value);
    }

    if (it != axis.begin()) {
        const std::size_t idx = static_cast<std::size_t>(std::distance(axis.begin(), it)) - 1;
        const double diff = std::abs(axis[idx] - value);
        if (diff < best_diff) {
            best = idx;
            best_diff = diff;
        }
    }

    if (!nearly_equal(axis[best], value, tol)) {
        return false;
    }

    out = static_cast<long long>(best);
    return true;
}

static bool extract_axis_aligned_cell(
    const Mesh& mesh,
    const HexElement& h,
    const ValidateOptions& opts,
    CellGeom& out,
    std::string& error
) {
    std::set<long long> distinct_nodes;
    for (long long nid : h.node_ids) {
        distinct_nodes.insert(nid);
    }

    if (distinct_nodes.size() != 8) {
        std::ostringstream oss;
        oss << "line " << h.line_number
            << ": element " << h.original_id
            << " does not contain 8 distinct node ids";
        error = oss.str();
        return false;
    }

    std::vector<Vec3> pts;
    pts.reserve(8);

    for (long long nid : h.node_ids) {
        const auto it = mesh.nodes.find(nid);
        if (it == mesh.nodes.end()) {
            std::ostringstream oss;
            oss << "line " << h.line_number
                << ": element " << h.original_id
                << " references missing node id " << nid;
            error = oss.str();
            return false;
        }
        pts.push_back(it->second.p);
    }

    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;
    xs.reserve(8);
    ys.reserve(8);
    zs.reserve(8);

    for (const Vec3& p : pts) {
        xs.push_back(p.x);
        ys.push_back(p.y);
        zs.push_back(p.z);
    }

    xs = unique_sorted_tol(xs, opts.tol);
    ys = unique_sorted_tol(ys, opts.tol);
    zs = unique_sorted_tol(zs, opts.tol);

    if (xs.size() != 2 || ys.size() != 2 || zs.size() != 2) {
        std::ostringstream oss;
        oss << "line " << h.line_number
            << ": element " << h.original_id
            << " is not an axis-aligned rectangular cell"
            << " because it has "
            << xs.size() << " unique x-values, "
            << ys.size() << " unique y-values, and "
            << zs.size() << " unique z-values";
        error = oss.str();
        return false;
    }

    const double xmin = xs[0];
    const double xmax = xs[1];
    const double ymin = ys[0];
    const double ymax = ys[1];
    const double zmin = zs[0];
    const double zmax = zs[1];

    const double dx = xmax - xmin;
    const double dy = ymax - ymin;
    const double dz = zmax - zmin;

    if (!positive_length(dx, opts.tol) ||
        !positive_length(dy, opts.tol) ||
        !positive_length(dz, opts.tol)) {
        std::ostringstream oss;
        oss << "line " << h.line_number
            << ": element " << h.original_id
            << " has non-positive or numerically degenerate extent";
        error = oss.str();
        return false;
    }

    bool seen[2][2][2] = {};
    for (const Vec3& p : pts) {
        const int ix = classify_to_pair(p.x, xmin, xmax, opts.tol);
        const int iy = classify_to_pair(p.y, ymin, ymax, opts.tol);
        const int iz = classify_to_pair(p.z, zmin, zmax, opts.tol);

        if (ix < 0 || iy < 0 || iz < 0) {
            std::ostringstream oss;
            oss << "line " << h.line_number
                << ": element " << h.original_id
                << " contains a node that is not located on one of the 8 axis-aligned rectangular corners";
            error = oss.str();
            return false;
        }

        seen[ix][iy][iz] = true;
    }

    for (int ix = 0; ix < 2; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int iz = 0; iz < 2; ++iz) {
                if (!seen[ix][iy][iz]) {
                    std::ostringstream oss;
                    oss << "line " << h.line_number
                        << ": element " << h.original_id
                        << " does not contain all 8 rectangular/cubical corner combinations";
                    error = oss.str();
                    return false;
                }
            }
        }
    }

    out.element_id = h.original_id;
    out.line_number = h.line_number;
    out.has_material = h.has_material;
    out.material = h.material;
    out.xmin = xmin;
    out.xmax = xmax;
    out.ymin = ymin;
    out.ymax = ymax;
    out.zmin = zmin;
    out.zmax = zmax;
    out.dx = dx;
    out.dy = dy;
    out.dz = dz;

    return true;
}

static void update_bounds(ValidationResult& result, const CellGeom& c) {
    if (!result.has_bounds) {
        result.min_bounds = Vec3{c.xmin, c.ymin, c.zmin};
        result.max_bounds = Vec3{c.xmax, c.ymax, c.zmax};
        result.has_bounds = true;
        return;
    }

    result.min_bounds.x = std::min(result.min_bounds.x, c.xmin);
    result.min_bounds.y = std::min(result.min_bounds.y, c.ymin);
    result.min_bounds.z = std::min(result.min_bounds.z, c.zmin);

    result.max_bounds.x = std::max(result.max_bounds.x, c.xmax);
    result.max_bounds.y = std::max(result.max_bounds.y, c.ymax);
    result.max_bounds.z = std::max(result.max_bounds.z, c.zmax);
}


static double percentile_sorted(const std::vector<double>& sorted, double p) {
    if (sorted.empty()) {
        return 0.0;
    }
    if (sorted.size() == 1U) {
        return sorted.front();
    }

    const double pos = p * static_cast<double>(sorted.size() - 1U);
    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = static_cast<std::size_t>(std::ceil(pos));
    const double frac = pos - static_cast<double>(lo);
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

static ScalarStats compute_scalar_stats(std::vector<double> values) {
    ScalarStats stats;
    stats.count = values.size();
    if (values.empty()) {
        return stats;
    }

    std::sort(values.begin(), values.end());
    stats.min = values.front();
    stats.max = values.back();
    stats.mean = std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
    stats.q1 = percentile_sorted(values, 0.25);
    stats.median = percentile_sorted(values, 0.50);
    stats.q3 = percentile_sorted(values, 0.75);

    double variance = 0.0;
    for (double v : values) {
        const double d = v - stats.mean;
        variance += d * d;
    }
    stats.stddev = std::sqrt(variance / static_cast<double>(values.size()));
    return stats;
}

static double cell_rel_spread(const CellGeom& c) {
    const double s_min = std::min(c.dx, std::min(c.dy, c.dz));
    const double s_max = std::max(c.dx, std::max(c.dy, c.dz));
    const double s_mean = (c.dx + c.dy + c.dz) / 3.0;
    return s_mean > 0.0 ? (s_max - s_min) / s_mean : std::numeric_limits<double>::infinity();
}

static SideLengthStats compute_side_length_stats(const std::vector<CellGeom>& cells) {
    SideLengthStats stats;
    if (cells.empty()) {
        return stats;
    }
    stats.available = true;

    std::vector<double> dxs;
    std::vector<double> dys;
    std::vector<double> dzs;
    std::vector<double> rels;
    dxs.reserve(cells.size());
    dys.reserve(cells.size());
    dzs.reserve(cells.size());
    rels.reserve(cells.size());

    double sum_abs_dx_dy = 0.0;
    double sum_abs_dx_dz = 0.0;
    double sum_abs_dy_dz = 0.0;
    std::vector<WorstCell> worst;
    worst.reserve(cells.size());

    for (const CellGeom& c : cells) {
        dxs.push_back(c.dx);
        dys.push_back(c.dy);
        dzs.push_back(c.dz);

        const double s_min = std::min(c.dx, std::min(c.dy, c.dz));
        const double s_max = std::max(c.dx, std::max(c.dy, c.dz));
        const double rel = cell_rel_spread(c);
        rels.push_back(rel);
        stats.max_abs_spread = std::max(stats.max_abs_spread, s_max - s_min);

        if (rel > 0.001) { ++stats.above_0_1_percent; }
        if (rel > 0.005) { ++stats.above_0_5_percent; }
        if (rel > 0.010) { ++stats.above_1_0_percent; }
        if (rel > 0.020) { ++stats.above_2_0_percent; }
        if (rel > 0.050) { ++stats.above_5_0_percent; }

        const double abs_dx_dy = std::abs(c.dx - c.dy);
        const double abs_dx_dz = std::abs(c.dx - c.dz);
        const double abs_dy_dz = std::abs(c.dy - c.dz);
        stats.max_abs_dx_dy = std::max(stats.max_abs_dx_dy, abs_dx_dy);
        stats.max_abs_dx_dz = std::max(stats.max_abs_dx_dz, abs_dx_dz);
        stats.max_abs_dy_dz = std::max(stats.max_abs_dy_dz, abs_dy_dz);
        sum_abs_dx_dy += abs_dx_dy;
        sum_abs_dx_dz += abs_dx_dz;
        sum_abs_dy_dz += abs_dy_dz;

        worst.push_back(WorstCell{c.element_id, c.line_number, c.dx, c.dy, c.dz, s_min, s_max, rel});
    }

    stats.dx = compute_scalar_stats(dxs);
    stats.dy = compute_scalar_stats(dys);
    stats.dz = compute_scalar_stats(dzs);
    stats.rel_spread = compute_scalar_stats(rels);
    stats.mean_abs_dx_dy = sum_abs_dx_dy / static_cast<double>(cells.size());
    stats.mean_abs_dx_dz = sum_abs_dx_dz / static_cast<double>(cells.size());
    stats.mean_abs_dy_dz = sum_abs_dy_dz / static_cast<double>(cells.size());

    std::sort(
        worst.begin(),
        worst.end(),
        [](const WorstCell& a, const WorstCell& b) {
            if (a.rel_spread != b.rel_spread) {
                return a.rel_spread > b.rel_spread;
            }
            return a.element_id < b.element_id;
        }
    );
    if (worst.size() > 10U) {
        worst.resize(10U);
    }
    stats.worst_cells = worst;
    return stats;
}

static AxisIntervalStats compute_axis_interval_stats(const std::vector<double>& axis, double tol) {
    AxisIntervalStats stats;
    stats.available = true;
    stats.coordinates = axis.size();
    if (axis.size() < 2U) {
        return stats;
    }

    std::vector<double> intervals;
    intervals.reserve(axis.size() - 1U);
    for (std::size_t i = 1; i < axis.size(); ++i) {
        intervals.push_back(axis[i] - axis[i - 1U]);
    }

    stats.intervals = intervals.size();
    const ScalarStats interval_stats = compute_scalar_stats(intervals);
    stats.min_interval = interval_stats.min;
    stats.max_interval = interval_stats.max;
    stats.mean_interval = interval_stats.mean;
    stats.median_interval = interval_stats.median;
    stats.unique_intervals = unique_sorted_tol(intervals, tol).size();
    return stats;
}

static void build_rectilinear_axes(
    const std::vector<CellGeom>& cells,
    double tol,
    std::vector<double>& x_axis,
    std::vector<double>& y_axis,
    std::vector<double>& z_axis
) {
    x_axis.clear();
    y_axis.clear();
    z_axis.clear();
    x_axis.reserve(cells.size() * 2U);
    y_axis.reserve(cells.size() * 2U);
    z_axis.reserve(cells.size() * 2U);

    for (const CellGeom& c : cells) {
        x_axis.push_back(c.xmin);
        x_axis.push_back(c.xmax);
        y_axis.push_back(c.ymin);
        y_axis.push_back(c.ymax);
        z_axis.push_back(c.zmin);
        z_axis.push_back(c.zmax);
    }

    x_axis = unique_sorted_tol(x_axis, tol);
    y_axis = unique_sorted_tol(y_axis, tol);
    z_axis = unique_sorted_tol(z_axis, tol);
}

static RectilinearAxisStats compute_rectilinear_axis_stats(
    const std::vector<double>& x_axis,
    const std::vector<double>& y_axis,
    const std::vector<double>& z_axis,
    double tol
) {
    RectilinearAxisStats stats;
    stats.available = true;
    stats.x = compute_axis_interval_stats(x_axis, tol);
    stats.y = compute_axis_interval_stats(y_axis, tol);
    stats.z = compute_axis_interval_stats(z_axis, tol);
    stats.median_interval = Vec3{stats.x.median_interval, stats.y.median_interval, stats.z.median_interval};
    return stats;
}

static bool map_rectilinear_cells(
    const std::vector<CellGeom>& cells,
    const std::vector<double>& x_axis,
    const std::vector<double>& y_axis,
    const std::vector<double>& z_axis,
    const ValidateOptions& opts,
    ValidationResult& result,
    std::unordered_set<CellIndex, CellIndexHash, CellIndexEqual>& occupied
) {
    bool ok_all = true;
    occupied.clear();
    occupied.reserve(cells.size() * 2U + 1U);

    for (const CellGeom& c : cells) {
        long long i0 = 0;
        long long i1 = 0;
        long long j0 = 0;
        long long j1 = 0;
        long long k0 = 0;
        long long k1 = 0;

        const bool ok =
            find_axis_index(x_axis, c.xmin, opts.tol, i0) &&
            find_axis_index(x_axis, c.xmax, opts.tol, i1) &&
            find_axis_index(y_axis, c.ymin, opts.tol, j0) &&
            find_axis_index(y_axis, c.ymax, opts.tol, j1) &&
            find_axis_index(z_axis, c.zmin, opts.tol, k0) &&
            find_axis_index(z_axis, c.zmax, opts.tol, k1);

        if (!ok) {
            std::ostringstream oss;
            oss << "line " << c.line_number
                << ": element " << c.element_id
                << " could not be mapped to rectilinear coordinate axes";
            add_validation_error(result, oss.str(), opts.max_errors);
            ok_all = false;
            continue;
        }

        if (i1 != i0 + 1 || j1 != j0 + 1 || k1 != k0 + 1) {
            std::ostringstream oss;
            oss << "line " << c.line_number
                << ": element " << c.element_id
                << " spans more than one rectilinear grid interval"
                << " and is therefore incompatible with a bitmap cubical complex";
            add_validation_error(result, oss.str(), opts.max_errors);
            ok_all = false;
            continue;
        }

        const CellIndex idx{i0, j0, k0};
        if (occupied.find(idx) != occupied.end()) {
            std::ostringstream oss;
            oss << "line " << c.line_number
                << ": element " << c.element_id
                << " duplicates occupied rectilinear cell " << cell_key_string(idx);
            add_validation_error(result, oss.str(), opts.max_errors);
            ok_all = false;
            continue;
        }

        occupied.insert(idx);
    }

    return ok_all;
}

static bool rectilinear_compatible_without_reporting(
    const std::vector<CellGeom>& cells,
    const ValidateOptions& opts
) {
    std::vector<double> x_axis;
    std::vector<double> y_axis;
    std::vector<double> z_axis;
    build_rectilinear_axes(cells, opts.tol, x_axis, y_axis, z_axis);

    ValidationResult scratch;
    std::unordered_set<CellIndex, CellIndexHash, CellIndexEqual> occupied;
    return map_rectilinear_cells(cells, x_axis, y_axis, z_axis, opts, scratch, occupied) &&
        occupied.size() == cells.size();
}

static void compute_face_statistics(
    const std::unordered_set<CellIndex, CellIndexHash, CellIndexEqual>& occupied,
    ValidationResult& result
) {
    std::size_t pairs = 0;

    for (const CellIndex& c : occupied) {
        const CellIndex xp{c.i + 1, c.j, c.k};
        const CellIndex yp{c.i, c.j + 1, c.k};
        const CellIndex zp{c.i, c.j, c.k + 1};

        if (occupied.find(xp) != occupied.end()) {
            ++pairs;
        }
        if (occupied.find(yp) != occupied.end()) {
            ++pairs;
        }
        if (occupied.find(zp) != occupied.end()) {
            ++pairs;
        }
    }

    result.face_adjacencies = pairs;
    result.boundary_faces = occupied.size() * 6U - 2U * pairs;
}

static ValidationResult validate_mesh(const Mesh& mesh, const ValidateOptions& opts) {
    ValidationResult result;

    result.n_nodes = mesh.nodes.size();
    result.n_hexes = mesh.hexes.size();
    result.n_blocks = mesh.blocks.size();

    for (const MeshBlock& b : mesh.blocks) {
        if (b.supported_hexa) {
            ++result.n_supported_blocks;
        } else {
            ++result.n_unsupported_blocks;
        }
    }

    result.warnings = mesh.warnings;

    if (opts.tol <= 0.0 || !std::isfinite(opts.tol)) {
        add_validation_error(result, "coordinate tolerance must be a finite positive number", opts.max_errors);
        result.valid = false;
        return result;
    }

    if (opts.cube_rel_tol < 0.0 || !std::isfinite(opts.cube_rel_tol)) {
        add_validation_error(result, "cube relative tolerance must be a finite non-negative number", opts.max_errors);
        result.valid = false;
        return result;
    }

    if (mesh.blocks.empty()) {
        add_validation_error(result, "no GiD mesh blocks were found", opts.max_errors);
    }

    if (mesh.hexes.empty()) {
        add_validation_error(result, "no supported 8-node hexahedral elements were found", opts.max_errors);
    }

    std::vector<CellGeom> cells;
    cells.reserve(mesh.hexes.size());

    for (const HexElement& h : mesh.hexes) {
        CellGeom c;
        std::string err;

        if (!extract_axis_aligned_cell(mesh, h, opts, c, err)) {
            add_validation_error(result, err, opts.max_errors);
            continue;
        }

        if (c.has_material) {
            result.material_counts[c.material] += 1;
        }

        update_bounds(result, c);
        cells.push_back(c);
    }

    result.side_stats = compute_side_length_stats(cells);

    if (!result.errors.empty()) {
        result.valid = false;
        return result;
    }

    if (cells.empty()) {
        add_validation_error(result, "no valid cells remain after geometric validation", opts.max_errors);
        result.valid = false;
        return result;
    }

    result.topology_compatible = rectilinear_compatible_without_reporting(cells, opts);

    std::unordered_set<CellIndex, CellIndexHash, CellIndexEqual> occupied;
    occupied.reserve(cells.size() * 2U + 1U);

    if (opts.grid_mode == GridMode::StrictCube) {
        const double h = cells.front().dx;
        std::size_t failing_cube_equality = 0;
        std::size_t failing_reference_spacing = 0;

        for (const CellGeom& c : cells) {
            if (!nearly_equal(c.dx, c.dy, opts.tol) ||
                !nearly_equal(c.dx, c.dz, opts.tol) ||
                !nearly_equal(c.dy, c.dz, opts.tol)) {
                ++failing_cube_equality;
            }

            if (!nearly_equal(c.dx, h, opts.tol) ||
                !nearly_equal(c.dy, h, opts.tol) ||
                !nearly_equal(c.dz, h, opts.tol)) {
                ++failing_reference_spacing;
            }
        }

        if (failing_cube_equality > 0U || failing_reference_spacing > 0U) {
            std::ostringstream oss;
            oss << failing_cube_equality
                << " cells fail strict cube edge equality and "
                << failing_reference_spacing
                << " cells differ from reference spacing h=" << h
                << "; maximum relative side-length spread="
                << result.side_stats.rel_spread.max
                << ", maximum absolute spread=" << result.side_stats.max_abs_spread
                << "; see cube-shape diagnostics and worst-cell table.";
            result.errors.push_back(oss.str());
            result.valid = false;
            return result;
        }

        result.origin = result.min_bounds;
        result.spacing = Vec3{h, h, h};

        long long max_i = 0;
        long long max_j = 0;
        long long max_k = 0;

        for (const CellGeom& c : cells) {
            long long i0 = 0;
            long long i1 = 0;
            long long j0 = 0;
            long long j1 = 0;
            long long k0 = 0;
            long long k1 = 0;

            const bool ok =
                coordinate_to_strict_index(c.xmin, result.origin.x, h, opts.tol, i0) &&
                coordinate_to_strict_index(c.xmax, result.origin.x, h, opts.tol, i1) &&
                coordinate_to_strict_index(c.ymin, result.origin.y, h, opts.tol, j0) &&
                coordinate_to_strict_index(c.ymax, result.origin.y, h, opts.tol, j1) &&
                coordinate_to_strict_index(c.zmin, result.origin.z, h, opts.tol, k0) &&
                coordinate_to_strict_index(c.zmax, result.origin.z, h, opts.tol, k1);

            if (!ok) {
                std::ostringstream oss;
                oss << "line " << c.line_number
                    << ": element " << c.element_id
                    << " is not aligned to the inferred strict lattice";
                add_validation_error(result, oss.str(), opts.max_errors);
                continue;
            }

            if (i1 != i0 + 1 || j1 != j0 + 1 || k1 != k0 + 1) {
                std::ostringstream oss;
                oss << "line " << c.line_number
                    << ": element " << c.element_id
                    << " does not occupy exactly one strict lattice cell";
                add_validation_error(result, oss.str(), opts.max_errors);
                continue;
            }

            const CellIndex idx{i0, j0, k0};

            if (occupied.find(idx) != occupied.end()) {
                std::ostringstream oss;
                oss << "line " << c.line_number
                    << ": element " << c.element_id
                    << " duplicates occupied lattice cell " << cell_key_string(idx);
                add_validation_error(result, oss.str(), opts.max_errors);
                continue;
            }

            occupied.insert(idx);

            max_i = std::max(max_i, i1);
            max_j = std::max(max_j, j1);
            max_k = std::max(max_k, k1);
        }

        result.nx = max_i;
        result.ny = max_j;
        result.nz = max_k;
    } else {
        std::vector<double> x_axis;
        std::vector<double> y_axis;
        std::vector<double> z_axis;
        build_rectilinear_axes(cells, opts.tol, x_axis, y_axis, z_axis);
        result.axis_stats = compute_rectilinear_axis_stats(x_axis, y_axis, z_axis, opts.tol);

        if (opts.grid_mode == GridMode::ApproximateCube) {
            std::size_t exceeding = 0;
            for (const CellGeom& c : cells) {
                if (cell_rel_spread(c) > opts.cube_rel_tol) {
                    ++exceeding;
                }
            }

            if (exceeding > 0U) {
                std::ostringstream oss;
                const double pct = 100.0 * static_cast<double>(exceeding) / static_cast<double>(cells.size());
                oss << exceeding << " cells exceed cube_rel_tol = " << opts.cube_rel_tol
                    << " (" << pct << "% of cells); see cube-shape diagnostics and worst-cell table.";
                result.errors.push_back(oss.str());
                result.valid = false;
                return result;
            }
        }

        result.origin = Vec3{x_axis.front(), y_axis.front(), z_axis.front()};
        result.spacing = Vec3{
            x_axis.size() > 1U ? x_axis[1] - x_axis[0] : 0.0,
            y_axis.size() > 1U ? y_axis[1] - y_axis[0] : 0.0,
            z_axis.size() > 1U ? z_axis[1] - z_axis[0] : 0.0
        };

        result.nx = static_cast<long long>(x_axis.size()) - 1;
        result.ny = static_cast<long long>(y_axis.size()) - 1;
        result.nz = static_cast<long long>(z_axis.size()) - 1;

        map_rectilinear_cells(cells, x_axis, y_axis, z_axis, opts, result, occupied);
    }

    if (!result.errors.empty()) {
        result.valid = false;
        return result;
    }

    result.occupied_cubes = occupied.size();
    compute_face_statistics(occupied, result);

    if (result.occupied_cubes != mesh.hexes.size()) {
        std::ostringstream oss;
        oss << "occupied cell count (" << result.occupied_cubes
            << ") differs from parsed hexahedral element count (" << mesh.hexes.size() << ")";
        add_validation_error(result, oss.str(), opts.max_errors);
    }

    result.valid = result.errors.empty();
    return result;
}

static std::string grid_mode_name(GridMode mode) {
    switch (mode) {
        case GridMode::StrictCube:
            return "strict-cube";
        case GridMode::ApproximateCube:
            return "approximate-cube";
        case GridMode::Rectilinear:
            return "rectilinear";
    }
    return "unknown";
}

static void print_version() {
    std::cout << "captop " << CAPTOP_VERSION << "\n";
    std::cout << "C++ standard target: C++17\n";
    std::cout << "GUDHI integration: not enabled in stages 1-3\n";
}

static void print_usage(std::ostream& os) {
    os << "CAPTOP - Cubical Analysis Pipeline for Topology\n\n";
    os << "Usage:\n";
    os << "  captop --version\n";
    os << "  captop validate <input.msh> [options]\n\n";
    os << "Commands:\n";
    os << "  validate                 Parse GiD ASCII mesh and validate cubic-grid compatibility\n\n";
    os << "Options for validate:\n";
    os << "  --grid strict-cube       Exact equal-edge cubes on a uniform lattice [default]\n";
    os << "  --grid approximate-cube  Axis-aligned cells nearly cubic within --cube-rel-tol\n";
    os << "  --grid rectilinear       Axis-aligned rectangular cells on rectilinear axes for topology-first analysis\n";
    os << "                            (legacy alias: --grid strict maps to strict-cube)\n";
    os << "  --tol <value>            Coordinate/snapping tolerance [default: 1e-8]\n";
    os << "  --cube-rel-tol <value>   Relative side-length spread tolerated in approximate-cube mode [default: 0.02]\n";
    os << "  --ignore-non-hexa        Ignore unsupported non-hexahedral mesh blocks\n";
    os << "  --max-errors <N>         Maximum number of validation errors to print [default: 30]\n";
    os << "  -h, --help               Show this help message\n";
}


static void print_scalar_stats_line(const std::string& label, const ScalarStats& s) {
    std::cout << "  " << label << " min/Q1/median/mean/Q3/max/std : "
              << s.min << " / " << s.q1 << " / " << s.median << " / "
              << s.mean << " / " << s.q3 << " / " << s.max << " / "
              << s.stddev << "\n";
}

static void print_axis_interval_line(const std::string& label, const AxisIntervalStats& s) {
    std::cout << "  " << label
              << " coordinates/intervals/min/median/max/mean/unique : "
              << s.coordinates << " / " << s.intervals << " / "
              << s.min_interval << " / " << s.median_interval << " / "
              << s.max_interval << " / " << s.mean_interval << " / "
              << s.unique_intervals << "\n";
}

static void print_validation_report(
    const std::string& input_path,
    const ValidateOptions& opts,
    const ValidationResult& r
) {
    std::cout << "============================================================\n";
    std::cout << "CAPTOP validation report\n";
    std::cout << "============================================================\n";
    std::cout << "Software version      : " << CAPTOP_VERSION << "\n";
    std::cout << "Input file            : " << input_path << "\n";
    std::cout << "Grid mode             : " << grid_mode_name(opts.grid_mode) << "\n";
    std::cout << "Coordinate tolerance  : " << std::setprecision(12) << opts.tol << "\n";
    if (opts.grid_mode == GridMode::ApproximateCube) {
        std::cout << "Cube relative tol     : " << opts.cube_rel_tol << "\n";
    }
    std::cout << "Status under selected mode: " << (r.valid ? "VALID" : "INVALID") << "\n";
    std::cout << "\n";

    std::cout << "GiD mesh structure\n";
    std::cout << "  Mesh blocks         : " << r.n_blocks << "\n";
    std::cout << "  Supported hex blocks: " << r.n_supported_blocks << "\n";
    std::cout << "  Unsupported blocks  : " << r.n_unsupported_blocks << "\n";
    std::cout << "  Nodes               : " << r.n_nodes << "\n";
    std::cout << "  Hexahedra parsed    : " << r.n_hexes << "\n";
    std::cout << "\n";

    if (r.has_bounds) {
        std::cout << "Coordinate bounds\n";
        std::cout << "  min                 : "
                  << r.min_bounds.x << " "
                  << r.min_bounds.y << " "
                  << r.min_bounds.z << "\n";
        std::cout << "  max                 : "
                  << r.max_bounds.x << " "
                  << r.max_bounds.y << " "
                  << r.max_bounds.z << "\n";
        std::cout << "\n";
    }

    if (r.side_stats.available) {
        std::cout << "Element side-length statistics\n";
        print_scalar_stats_line("dx", r.side_stats.dx);
        print_scalar_stats_line("dy", r.side_stats.dy);
        print_scalar_stats_line("dz", r.side_stats.dz);
        std::cout << "\n";

        std::cout << "Approximate cube diagnostics\n";
        std::cout << "  relative spread min/median/mean/max : "
                  << r.side_stats.rel_spread.min << " / "
                  << r.side_stats.rel_spread.median << " / "
                  << r.side_stats.rel_spread.mean << " / "
                  << r.side_stats.rel_spread.max << "\n";
        std::cout << "  relative spread Q1/Q3              : "
                  << r.side_stats.rel_spread.q1 << " / "
                  << r.side_stats.rel_spread.q3 << "\n";
        std::cout << "  cells above 0.1% / 0.5% / 1% / 2% / 5% : "
                  << r.side_stats.above_0_1_percent << " / "
                  << r.side_stats.above_0_5_percent << " / "
                  << r.side_stats.above_1_0_percent << " / "
                  << r.side_stats.above_2_0_percent << " / "
                  << r.side_stats.above_5_0_percent << "\n";
        std::cout << "  max abs spread                    : " << r.side_stats.max_abs_spread << "\n";
        std::cout << "  mean abs(dx-dy/dx-dz/dy-dz)       : "
                  << r.side_stats.mean_abs_dx_dy << " / "
                  << r.side_stats.mean_abs_dx_dz << " / "
                  << r.side_stats.mean_abs_dy_dz << "\n";
        std::cout << "  max abs(dx-dy/dx-dz/dy-dz)        : "
                  << r.side_stats.max_abs_dx_dy << " / "
                  << r.side_stats.max_abs_dx_dz << " / "
                  << r.side_stats.max_abs_dy_dz << "\n";
        std::cout << "\n";

        if (!r.side_stats.worst_cells.empty()) {
            std::cout << "Worst cells by side-length relative spread\n";
            std::cout << "  line element dx dy dz s_min s_max rel_spread_percent\n";
            for (const WorstCell& w : r.side_stats.worst_cells) {
                std::cout << "  " << w.line_number << " " << w.element_id << " "
                          << w.dx << " " << w.dy << " " << w.dz << " "
                          << w.s_min << " " << w.s_max << " "
                          << (100.0 * w.rel_spread) << "\n";
            }
            std::cout << "\n";
        }
    }

    if (r.axis_stats.available) {
        std::cout << "Rectilinear axis interval statistics\n";
        print_axis_interval_line("x", r.axis_stats.x);
        print_axis_interval_line("y", r.axis_stats.y);
        print_axis_interval_line("z", r.axis_stats.z);
        std::cout << "\n";
    }

    if (r.valid) {
        std::cout << "Cubical grid\n";
        std::cout << "  Grid dimensions     : "
                  << r.nx << " x " << r.ny << " x " << r.nz << "\n";
        std::cout << "  Origin              : "
                  << r.origin.x << " "
                  << r.origin.y << " "
                  << r.origin.z << "\n";

        if (opts.grid_mode == GridMode::StrictCube) {
            std::cout << "  Spacing             : "
                      << r.spacing.x << " "
                      << r.spacing.y << " "
                      << r.spacing.z << "\n";
        } else if (opts.grid_mode == GridMode::ApproximateCube) {
            std::cout << "  Representative spacing : "
                      << r.side_stats.dx.median << " "
                      << r.side_stats.dy.median << " "
                      << r.side_stats.dz.median << "\n";
            std::cout << "  First axis interval : "
                      << r.spacing.x << " "
                      << r.spacing.y << " "
                      << r.spacing.z << "\n";
            std::cout << "  Median axis interval: "
                      << r.axis_stats.median_interval.x << " "
                      << r.axis_stats.median_interval.y << " "
                      << r.axis_stats.median_interval.z << "\n";
        } else {
            std::cout << "  First axis interval : "
                      << r.spacing.x << " "
                      << r.spacing.y << " "
                      << r.spacing.z << "\n";
            std::cout << "  Median axis interval: "
                      << r.axis_stats.median_interval.x << " "
                      << r.axis_stats.median_interval.y << " "
                      << r.axis_stats.median_interval.z << "\n";
        }

        std::cout << "  Occupied cells      : " << r.occupied_cubes << "\n";
        std::cout << "  Face adjacencies    : " << r.face_adjacencies << "\n";
        std::cout << "  Boundary faces      : " << r.boundary_faces << "\n";
        std::cout << "\n";
    }

    if (!r.material_counts.empty()) {
        std::cout << "Materials/layers\n";
        for (const auto& kv : r.material_counts) {
            std::cout << "  " << kv.first << " : " << kv.second << " elements\n";
        }
        std::cout << "\n";
    }

    if (!r.warnings.empty()) {
        std::cout << "Warnings\n";
        for (const std::string& w : r.warnings) {
            std::cout << "  - " << w << "\n";
        }
        std::cout << "\n";
    }

    if (!r.errors.empty()) {
        std::cout << "Errors\n";
        for (const std::string& e : r.errors) {
            std::cout << "  - " << e << "\n";
        }
        std::cout << "\n";
    }

    std::cout << "Scientific interpretation\n";
    if (r.valid) {
        if (opts.grid_mode == GridMode::StrictCube) {
            std::cout << "  Exact metric conclusion: exact-cube compatible under the selected tolerance.\n";
        } else if (opts.grid_mode == GridMode::ApproximateCube) {
            std::cout << "  Approximate geometric conclusion: approximately cubic within the selected relative tolerance.\n";
        } else {
            std::cout << "  Topological conclusion: rectilinear-compatible with well-defined face adjacencies.\n";
        }
    } else if (r.topology_compatible) {
        if (opts.grid_mode == GridMode::StrictCube) {
            std::cout << "  The mesh is axis-aligned and rectilinear-compatible, but not exact-cube compatible under the selected tolerance.\n";
            std::cout << "  Approximate geometric status may be valid: rerun with --grid approximate-cube to apply --cube-rel-tol.\n";
        } else if (opts.grid_mode == GridMode::ApproximateCube) {
            std::cout << "  The mesh is axis-aligned and rectilinear-compatible, but not approximate-cube compatible under the selected relative tolerance.\n";
        } else {
            std::cout << "  The mesh is axis-aligned but failed the selected rectilinear policy; see errors.\n";
        }
        std::cout << "  Topology-compatible status: likely valid under rectilinear mode.\n";
        std::cout << "  Topology-compatible rectilinear grid likely: rerun with --grid rectilinear if topology-first analysis is intended.\n";
    } else {
        std::cout << "  The mesh does not currently validate under the selected mode; see errors for geometric or grid-mapping failures.\n";
    }
    std::cout << "\n";

    std::cout << "Ready for Stage 4 bitmap conversion under selected mode: "
              << (r.valid ? "yes" : "no") << "\n";
    std::cout << "============================================================\n";
}

static int run_validate(int argc, char** argv) {
    if (argc < 3) {
        print_usage(std::cerr);
        return 1;
    }

    if (std::string(argv[2]) == "--help" || std::string(argv[2]) == "-h") {
        print_usage(std::cout);
        return 0;
    }

    std::string input_path = argv[2];

    ParseOptions parse_opts;
    ValidateOptions validate_opts;

    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];

        if (arg == "--grid") {
            if (i + 1 >= argc) {
                std::cerr << "error: --grid requires a value: strict-cube, approximate-cube, or rectilinear\n";
                return 1;
            }

            const std::string value = to_lower(argv[++i]);
            if (value == "strict" || value == "strict-cube") {
                validate_opts.grid_mode = GridMode::StrictCube;
            } else if (value == "approximate-cube") {
                validate_opts.grid_mode = GridMode::ApproximateCube;
            } else if (value == "rectilinear") {
                validate_opts.grid_mode = GridMode::Rectilinear;
            } else {
                std::cerr << "error: unsupported grid mode '" << value << "'\n";
                return 1;
            }
        } else if (arg == "--tol") {
            if (i + 1 >= argc) {
                std::cerr << "error: --tol requires a numeric value\n";
                return 1;
            }

            double tol = 0.0;
            if (!parse_double(argv[++i], tol) || tol <= 0.0) {
                std::cerr << "error: invalid tolerance value\n";
                return 1;
            }

            validate_opts.tol = tol;
        } else if (arg == "--cube-rel-tol") {
            if (i + 1 >= argc) {
                std::cerr << "error: --cube-rel-tol requires a numeric value\n";
                return 1;
            }

            double cube_rel_tol = 0.0;
            if (!parse_double(argv[++i], cube_rel_tol) || cube_rel_tol < 0.0) {
                std::cerr << "error: invalid --cube-rel-tol value\n";
                return 1;
            }

            validate_opts.cube_rel_tol = cube_rel_tol;
        } else if (arg == "--ignore-non-hexa") {
            parse_opts.ignore_non_hexa = true;
        } else if (arg == "--max-errors") {
            if (i + 1 >= argc) {
                std::cerr << "error: --max-errors requires an integer value\n";
                return 1;
            }

            long long n = 0;
            if (!parse_long_long(argv[++i], n) || n <= 0) {
                std::cerr << "error: invalid --max-errors value\n";
                return 1;
            }

            validate_opts.max_errors = static_cast<std::size_t>(n);
        } else if (arg == "-h" || arg == "--help") {
            print_usage(std::cout);
            return 0;
        } else {
            std::cerr << "error: unknown option '" << arg << "'\n";
            return 1;
        }
    }

    try {
        Mesh mesh = parse_gid_mesh(input_path, parse_opts);
        ValidationResult result = validate_mesh(mesh, validate_opts);
        print_validation_report(input_path, validate_opts, result);
        return result.valid ? 0 : 2;
    } catch (const ParseError& e) {
        std::cerr << "parse error: " << e.what() << "\n";
        return 1;
    } catch (const ValidationError& e) {
        std::cerr << "validation error: " << e.what() << "\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "fatal error: " << e.what() << "\n";
        return 1;
    }
}

} // namespace captop

int main(int argc, char** argv) {
    if (argc <= 1) {
        captop::print_usage(std::cerr);
        return 1;
    }

    const std::string command = argv[1];

    if (command == "--version" || command == "-v") {
        captop::print_version();
        return 0;
    }

    if (command == "--help" || command == "-h") {
        captop::print_usage(std::cout);
        return 0;
    }

    if (command == "validate") {
        return captop::run_validate(argc, argv);
    }

    std::cerr << "error: unknown command '" << command << "'\n";
    captop::print_usage(std::cerr);
    return 1;
}

