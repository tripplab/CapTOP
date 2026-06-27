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
#include <chrono>
#include <exception>
#include <cstdint>
#include <filesystem>
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
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <queue>
#include <vector>
#ifdef CAPTOP_WITH_GUDHI
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>
#endif

namespace captop {

static const char* CAPTOP_VERSION = "0.1.0-stage6";

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

struct IndexedCell {
    long long i = 0;
    long long j = 0;
    long long k = 0;
    long long original_element_id = 0;
    bool has_material = false;
    long long material = 0;
    std::size_t line_number = 0;
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
    std::vector<IndexedCell> indexed_cells;
    std::vector<double> x_axis;
    std::vector<double> y_axis;
    std::vector<double> z_axis;
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


static std::string strip_octreemesh_comment(const std::string& line) {
    const std::size_t pos = line.find(';');
    if (pos == std::string::npos) {
        return trim(line);
    }
    return trim(line.substr(0, pos));
}

static std::vector<std::string> tokenize_plain_line(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string tok;
    while (iss >> tok) {
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


static std::string first_significant_line(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) {
        throw ParseError("cannot open input file '" + path + "'");
    }

    std::string line;
    while (std::getline(in, line)) {
        std::string stripped = trim(line);
        if (stripped.empty() || starts_with_case_insensitive(stripped, "#")) {
            continue;
        }
        const std::string without_octree_comment = strip_octreemesh_comment(stripped);
        if (without_octree_comment.empty()) {
            continue;
        }
        return without_octree_comment;
    }

    throw ParseError("input file '" + path + "' is empty or contains only comments");
}

static Mesh parse_octreemesh_file(const std::string& path, const ParseOptions&) {
    std::ifstream in(path.c_str());
    if (!in) {
        throw ParseError("cannot open input file '" + path + "'");
    }

    Mesh mesh;
    MeshBlock block;
    block.id = 0;
    block.dimension = 3;
    block.elemtype = "Hexahedra";
    block.nnode = 8;
    block.supported_hexa = true;
    mesh.blocks.push_back(block);

    auto next_data_line = [&](std::string& out, std::size_t& out_line) -> bool {
        std::string line;
        while (std::getline(in, line)) {
            ++out_line;
            const std::string stripped = strip_octreemesh_comment(line);
            if (stripped.empty()) {
                continue;
            }
            out = stripped;
            return true;
        }
        return false;
    };

    std::string line;
    std::size_t line_number = 0;

    if (!next_data_line(line, line_number) || line != "{Nodes}") {
        throw ParseError(line_context(line_number, "OctreeMesh file must start with {Nodes}"));
    }
    mesh.blocks[0].line_number = line_number;

    if (!next_data_line(line, line_number)) {
        throw ParseError("unexpected end of OctreeMesh file while reading node dimension");
    }
    int dimension = 0;
    {
        const std::vector<std::string> tokens = tokenize_plain_line(line);
        if (tokens.size() != 1 || !parse_int(tokens[0], dimension)) {
            throw ParseError(line_context(line_number, "OctreeMesh node dimension must be a single integer"));
        }
    }
    if (dimension != 3) {
        throw ParseError(line_context(line_number, "OctreeMesh dimension " + std::to_string(dimension) + " is not supported; only 3 is supported"));
    }
    mesh.blocks[0].dimension = dimension;

    if (!next_data_line(line, line_number)) {
        throw ParseError("unexpected end of OctreeMesh file while reading nodes count");
    }
    long long node_count_ll = 0;
    {
        const std::vector<std::string> tokens = tokenize_plain_line(line);
        if (tokens.size() != 1 || !parse_long_long(tokens[0], node_count_ll) || node_count_ll < 0) {
            throw ParseError(line_context(line_number, "OctreeMesh nodes count must be a non-negative integer"));
        }
    }
    const std::size_t node_count = static_cast<std::size_t>(node_count_ll);
    std::cerr << "Parsing OctreeMesh nodes: " << node_count << "\n";

    for (std::size_t i = 0; i < node_count; ++i) {
        if (!next_data_line(line, line_number)) {
            throw ParseError("unexpected end of OctreeMesh file while reading node " + std::to_string(i + 1));
        }
        const std::vector<std::string> tokens = tokenize_plain_line(line);
        if (tokens.size() != 3) {
            throw ParseError(line_context(line_number, "OctreeMesh node " + std::to_string(i + 1) + " must contain exactly 3 coordinates"));
        }
        double x = 0.0, y = 0.0, z = 0.0;
        if (!parse_double(tokens[0], x) || !parse_double(tokens[1], y) || !parse_double(tokens[2], z)) {
            throw ParseError(line_context(line_number, "OctreeMesh node " + std::to_string(i + 1) + " has invalid coordinate value"));
        }
        const long long node_id = static_cast<long long>(i + 1);
        Node n;
        n.original_id = node_id;
        n.p = Vec3{x, y, z};
        mesh.nodes.insert(std::make_pair(node_id, n));
    }
    std::cerr << "Finished OctreeMesh nodes: " << mesh.nodes.size() << "\n";

    if (!next_data_line(line, line_number) || line != "{Mesh}") {
        throw ParseError(line_context(line_number, "OctreeMesh {Mesh} section must follow the {Nodes} section"));
    }

    if (!next_data_line(line, line_number)) {
        throw ParseError("unexpected end of OctreeMesh file while reading element type");
    }
    int element_type = 0;
    {
        const std::vector<std::string> tokens = tokenize_plain_line(line);
        if (tokens.size() != 1 || !parse_int(tokens[0], element_type)) {
            throw ParseError(line_context(line_number, "OctreeMesh element type must be a single integer"));
        }
    }
    if (element_type != 5) {
        throw ParseError(line_context(line_number, "OctreeMesh element type " + std::to_string(element_type) + " is not supported; only 5=Hexahedra is supported"));
    }

    if (!next_data_line(line, line_number)) {
        throw ParseError("unexpected end of OctreeMesh file while reading nodes per element");
    }
    int nodes_per_element = 0;
    {
        const std::vector<std::string> tokens = tokenize_plain_line(line);
        if (tokens.size() != 1 || !parse_int(tokens[0], nodes_per_element)) {
            throw ParseError(line_context(line_number, "OctreeMesh nodes per element must be a single integer"));
        }
    }
    if (nodes_per_element != 8) {
        throw ParseError(line_context(line_number, "OctreeMesh nodes per element " + std::to_string(nodes_per_element) + " is not supported; only 8 is supported"));
    }
    mesh.blocks[0].nnode = nodes_per_element;

    if (!next_data_line(line, line_number)) {
        throw ParseError("unexpected end of OctreeMesh file while reading elements count");
    }
    long long element_count_ll = 0;
    {
        const std::vector<std::string> tokens = tokenize_plain_line(line);
        if (tokens.size() != 1 || !parse_long_long(tokens[0], element_count_ll) || element_count_ll < 0) {
            throw ParseError(line_context(line_number, "OctreeMesh elements count must be a non-negative integer"));
        }
    }
    const std::size_t element_count = static_cast<std::size_t>(element_count_ll);
    mesh.blocks[0].element_count = element_count;
    std::cerr << "Parsing OctreeMesh elements: " << element_count << "\n";

    for (std::size_t i = 0; i < element_count; ++i) {
        if (!next_data_line(line, line_number)) {
            throw ParseError("unexpected end of OctreeMesh file while reading element " + std::to_string(i + 1));
        }
        const std::vector<std::string> tokens = tokenize_plain_line(line);
        if (tokens.size() != 9) {
            throw ParseError(line_context(line_number, "OctreeMesh element " + std::to_string(i + 1) + " must contain exactly material plus 8 node ids"));
        }

        HexElement h;
        h.original_id = static_cast<long long>(i + 1);
        h.mesh_block_id = 0;
        h.line_number = line_number;
        h.has_material = true;
        if (!parse_long_long(tokens[0], h.material)) {
            throw ParseError(line_context(line_number, "OctreeMesh element " + std::to_string(i + 1) + " has invalid material id"));
        }
        for (std::size_t k = 0; k < 8; ++k) {
            long long nid = 0;
            if (!parse_long_long(tokens[1 + k], nid)) {
                throw ParseError(line_context(line_number, "OctreeMesh element " + std::to_string(i + 1) + " has invalid node id '" + tokens[1 + k] + "'"));
            }
            if (nid < 1 || nid > node_count_ll) {
                throw ParseError(line_context(line_number, "OctreeMesh element " + std::to_string(i + 1) + " references node id " + std::to_string(nid) + " outside valid range 1.." + std::to_string(node_count_ll)));
            }
            h.node_ids[k] = nid;
        }
        mesh.hexes.push_back(h);
    }
    std::cerr << "Finished OctreeMesh elements: " << mesh.hexes.size() << "\n";

    std::string extra;
    std::size_t extra_line = line_number;
    if (next_data_line(extra, extra_line)) {
        throw ParseError(line_context(extra_line, "unexpected data after OctreeMesh elements section"));
    }

    return mesh;
}

static Mesh parse_mesh_file(const std::string& path, const ParseOptions& opts) {
    const std::string first = first_significant_line(path);
    const std::string lower = to_lower(first);
    if (starts_with_case_insensitive(lower, "mesh")) {
        return parse_gid_mesh(path, opts);
    }
    if (first == "{Nodes}") {
        return parse_octreemesh_file(path, opts);
    }
    throw ParseError("unsupported mesh format in '" + path + "'; expected GiD 'mesh' header or OctreeMesh {Nodes} header");
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
    std::unordered_set<CellIndex, CellIndexHash, CellIndexEqual>& occupied,
    std::vector<IndexedCell>* indexed_cells = nullptr
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
        if (indexed_cells != nullptr) {
            indexed_cells->push_back(IndexedCell{i0, j0, k0, c.element_id, c.has_material, c.material, c.line_number});
        }
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
            result.indexed_cells.push_back(IndexedCell{i0, j0, k0, c.element_id, c.has_material, c.material, c.line_number});

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

        result.x_axis = x_axis;
        result.y_axis = y_axis;
        result.z_axis = z_axis;

        result.origin = Vec3{x_axis.front(), y_axis.front(), z_axis.front()};
        result.spacing = Vec3{
            x_axis.size() > 1U ? x_axis[1] - x_axis[0] : 0.0,
            y_axis.size() > 1U ? y_axis[1] - y_axis[0] : 0.0,
            z_axis.size() > 1U ? z_axis[1] - z_axis[0] : 0.0
        };

        result.nx = static_cast<long long>(x_axis.size()) - 1;
        result.ny = static_cast<long long>(y_axis.size()) - 1;
        result.nz = static_cast<long long>(z_axis.size()) - 1;

        map_rectilinear_cells(cells, x_axis, y_axis, z_axis, opts, result, occupied, &result.indexed_cells);
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
            return "strict";
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
#ifdef CAPTOP_WITH_GUDHI
    std::cout << "GUDHI support: enabled\n";
#else
    std::cout << "GUDHI support: disabled\n";
#endif
    std::cout << "Stage 5 fallback diagnostics: enabled\n";
}

static void print_validation_options(std::ostream& os) {
    os << "Validation options:\n";
    os << "  --grid strict-cube       Exact equal-edge cubes on a uniform lattice [default]\n";
    os << "  --grid approximate-cube  Axis-aligned cells nearly cubic within --cube-rel-tol\n";
    os << "  --grid rectilinear       Axis-aligned rectangular cells on rectilinear axes for topology-first analysis\n";
    os << "                            (legacy alias: --grid strict maps to strict-cube)\n";
    os << "  --tol <value>            Coordinate/snapping tolerance [default: 1e-8]\n";
    os << "  --cube-rel-tol <value>   Relative side-length spread tolerated in approximate-cube mode [default: 0.02]\n";
    os << "  --ignore-non-hexa        Ignore unsupported non-hexahedral mesh blocks\n";
    os << "  --max-errors <N>         Maximum number of validation errors to print [default: 30]\n";
}

static void print_convert_options(std::ostream& os) {
    os << "Convert options:\n";
    os << "  --out <dir>              Output directory [default: captop_out]\n";
    os << "  --filtration <policy>    Filtration policy: occupancy, material, scalar-file, or binary-threshold [default: occupancy]\n";
    os << "  --selected-material <N>  For material filtration, include only the selected material as finite\n";
    os << "  --scalar-file <csv>      CSV containing element_id,value columns for scalar-file or binary-threshold filtration\n";
    os << "  --threshold <value>      Threshold value required by binary-threshold filtration\n";
    os << "  --threshold-op <op>      Threshold operator: lt, le, gt, ge, eq, or ne [default: ge]\n";
    os << "  --missing-value inf      Missing cell value; only inf is supported\n";
    os << "  --overwrite              Replace existing CAPTOP output files in --out\n";
    os << "  --force                  Continue when the dense memory estimate exceeds --max-memory-gb\n";
    os << "  --max-memory-gb <value>  Dense-grid memory limit before --force is required [default: 4]\n";
}

static void print_persist_options(std::ostream& os) {
    os << "Persist options:\n";
    os << "  --homology-dim <dims>    Comma-separated dimensions from 0,1,2,3 [default: 0,1,2]\n";
    os << "  --field <prime>          Coefficient field [default: 2]\n";
    os << "  --filtration <policy>    Filtration policy: occupancy, material, or scalar-file [default: occupancy]\n";
    os << "  --scalar-file <csv>      CSV containing element_id,value columns for scalar-file filtration\n";
    os << "  --mode <mode>            Filtration mode: sublevel or superlevel [default: sublevel]\n";
    os << "  --min-persistence <eps>  Minimum persistence threshold [default: 0]\n";
    os << "  --out <dir>              Output directory [default: captop_persist_out]\n";
    os << "  --overwrite              Replace existing persistence output files\n";
    os << "  --force                  Continue when the dense memory estimate exceeds --max-memory-gb\n";
    os << "  --max-memory-gb <value>  Dense-grid memory limit before --force is required [default: 4]\n";
    os << "  --write-pairs            Write persistence_pairs.csv output\n";
    os << "  --write-diagrams         Write per-dimension diagram CSV outputs\n";
    os << "  --write-betti-curve      Write betti_curve.csv output\n";
    os << "  --write-barcode-summary  Write barcode_summary.csv output\n";
    os << "  --write-json             Write JSON summary output\n";
    os << "  --write-all              Enable all persistence output writers\n";
    os << "  --betti-curve-samples <N>  Number of uniform Betti-curve samples [default: 200]\n";
    os << "  --betti-curve-values <mode>  Betti-curve values: unique or uniform [default: unique]\n";
    os << "  --finite-only            Exclude essential intervals from persistence outputs\n";
    os << "  --include-essential      Include essential intervals in persistence outputs [default]\n";
    os << "  --quiet                  Suppress terminal report\n";
}

static void print_validate_usage(std::ostream& os) {
    os << "Usage:\n";
    os << "  captop validate <input-mesh> [options]\n\n";
    print_validation_options(os);
    os << "  -h, --help               Show this help message\n";
}

static void print_convert_usage(std::ostream& os) {
    os << "Usage:\n";
    os << "  captop convert <input-mesh> [options]\n";
    os << "  captop betti <input-mesh> [options]\n";
    os << "  captop persist <input-mesh> [options]\n\n";
    print_validation_options(os);
    os << "\n";
    print_convert_options(os);
    os << "\nBetti options:\n";
    os << "  --field <prime>          Coefficient field [default: 2]\n";
    os << "  --out <dir>              Output directory [default: captop_betti_out]\n";
    os << "  --overwrite              Replace existing Betti output files\n";
    os << "  --max-memory-gb <value>  Dense-grid memory limit [default: 2]\n";
    os << "  --force                  Continue above memory limit\n";
    os << "  --write-diagram          Write raw persistence intervals for debugging\n";
    os << "  --write-json             Write JSON summary (included by default)\n";
    os << "  --write-csv              Write CSV summary (included by default)\n";
    os << "  --quiet                  Suppress terminal report\n";
    os << "  -h, --help               Show this help message\n";
    os << "\n";
    print_persist_options(os);
}

static void print_usage(std::ostream& os) {
    os << "CAPTOP - Cubical Analysis Pipeline for Topology\n\n";
    os << "Usage:\n";
    os << "  captop --version\n";
    os << "  captop validate <input-mesh> [options]\n";
    os << "  captop convert <input-mesh> [options]\n";
    os << "  captop betti <input-mesh> [options]\n";
    os << "  captop persist <input-mesh> [options]\n\n";
    os << "Commands:\n";
    os << "  validate                 Parse GiD ASCII or OctreeMesh input and validate cubic-grid compatibility\n";
    os << "  convert                  Convert validated mesh to dense indexed cubical bitmap files\n";
    os << "  betti                   Compute Stage 5 occupied-domain Betti descriptors\n  persist                 Compute Stage 6 persistent homology with GUDHI\n\n";
    print_validation_options(os);
    os << "\n";
    print_convert_options(os);
    os << "\nBetti options:\n";
    os << "  --field <prime>          Coefficient field [default: 2]\n";
    os << "  --out <dir>              Output directory [default: captop_betti_out]\n";
    os << "  --overwrite              Replace existing Betti output files\n";
    os << "  --max-memory-gb <value>  Dense-grid memory limit [default: 2]\n";
    os << "  --force                  Continue above memory limit\n";
    os << "  --write-diagram          Write raw persistence intervals for debugging\n";
    os << "  --write-json             Write JSON summary (included by default)\n";
    os << "  --write-csv              Write CSV summary (included by default)\n";
    os << "  --quiet                  Suppress terminal report\n";
    os << "  -h, --help               Show this help message\n";
    os << "\n";
    print_persist_options(os);
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
        print_validate_usage(std::cout);
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
        Mesh mesh = parse_mesh_file(input_path, parse_opts);
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


struct GridCell {
    bool occupied = false;
    double cube_value = std::numeric_limits<double>::infinity();
    long long material = 0;
    bool has_material = false;
    long long original_element_id = -1;
};

struct ConvertOptions {
    std::string out_dir = "captop_out";
    std::string filtration = "occupancy";
    std::string scalar_file;
    bool has_threshold = false;
    double threshold = 0.0;
    std::string threshold_op = "ge";
    bool has_selected_material = false;
    long long selected_material = 0;
    bool force = false;
    bool overwrite = false;
    double max_memory_gb = 4.0;
};

struct ConvertCounts { std::size_t occupied=0, missing=0, finite=0, selected=0; };

static std::size_t linear_index(long long i,long long j,long long k,long long nx,long long ny){
    return static_cast<std::size_t>(i + nx * (j + ny * k));
}

static std::string json_escape(const std::string& v){
    std::ostringstream o; for(char c: v){ if(c=='"'||c=='\\') o<<'\\'<<c; else if(c=='\n') o<<"\\n"; else o<<c; } return o.str();
}
static void write_axis_json(std::ostream& out,const std::vector<double>& a){
    out << "["; for(size_t i=0;i<a.size();++i){ if(i) out << ", "; out << std::setprecision(17) << a[i]; } out << "]";
}

static std::unordered_map<long long,double> read_scalar_file(const std::string& path){
    std::ifstream in(path.c_str()); if(!in) throw std::runtime_error("cannot open scalar file '"+path+"'");
    std::string line; if(!std::getline(in,line)) throw std::runtime_error("scalar file is empty");
    auto split=[](const std::string& x){ std::vector<std::string> r; std::string cur; std::istringstream ss(x); while(std::getline(ss,cur,',')) r.push_back(trim(cur)); return r; };
    auto h=split(line); int eid=-1,val=-1; for(size_t i=0;i<h.size();++i){ if(h[i]=="element_id") eid=i; if(h[i]=="value") val=i; }
    if(eid<0||val<0) throw std::runtime_error("scalar file header must contain element_id and value columns");
    if(h.size()>2) std::cerr << "warning: extra scalar CSV columns are ignored\n";
    std::unordered_map<long long,double> m; size_t ln=1;
    while(std::getline(in,line)){ ++ln; if(trim(line).empty()) continue; auto f=split(line); if(f.size()<=static_cast<size_t>(std::max(eid,val))) throw std::runtime_error("scalar file line "+std::to_string(ln)+" has missing values");
        long long id; double v; if(!parse_long_long(f[eid],id)) throw std::runtime_error("scalar file line "+std::to_string(ln)+" has unparsable element_id");
        if(!parse_double(f[val],v)) throw std::runtime_error("scalar file line "+std::to_string(ln)+" has nonfinite or unparsable value");
        if (m.count(id)) {
            throw std::runtime_error("duplicate scalar value for element " + std::to_string(id));
        }
        m[id] = v;
    }
    return m;
}

static bool threshold_pass(double v,double t,const std::string& op,double tol){
    if (op == "lt") { return v < t; }
    if (op == "le") { return v <= t; }
    if (op == "gt") { return v > t; }
    if (op == "ge") { return v >= t; }
    if (op == "eq") { return nearly_equal(v, t, tol); }
    if (op == "ne") { return !nearly_equal(v, t, tol); }
    throw std::runtime_error("unknown threshold operator '" + op + "'");
}

template<class T> static void write_raw(const std::filesystem::path& path,const std::vector<T>& v){
    std::ofstream out(path, std::ios::binary); if(!out) throw std::runtime_error("cannot open output file '"+path.string()+"'");
    out.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size()*sizeof(T))); if(!out) throw std::runtime_error("could not write complete raw file '"+path.string()+"'");
}

static std::string conversion_report(const std::string& input,const std::string& outdir,const ValidateOptions& vopts,const ConvertOptions& copts,const ValidationResult& vr,const ConvertCounts& c,size_t bytes){
    std::ostringstream r; r<<"============================================================\nCAPTOP conversion report\n============================================================\n";
    r<<"Software version      : "<<CAPTOP_VERSION<<"\nInput file            : "<<input<<"\nOutput directory      : "<<outdir<<"\nGrid mode             : "<<grid_mode_name(vopts.grid_mode)<<"\nTolerance             : "<<std::setprecision(12)<<vopts.tol<<"\nFiltration policy     : "<<copts.filtration<<"\nStatus                : VALID AND CONVERTED\n\nGrid\n  Dimensions          : "<<vr.nx<<" x "<<vr.ny<<" x "<<vr.nz<<"\n  Total cells         : "<<(c.occupied+c.missing)<<"\n  Occupied cells      : "<<c.occupied<<"\n  Missing cells       : "<<c.missing<<"\n  Finite-value cells  : "<<c.finite<<"\n  Selected cells      : "<<c.selected<<"\n  Dense memory estimate: "<<bytes<<" bytes\n\nFiles written\n  Metadata            : "<<outdir<<"/captop_grid_metadata.json\n  Cube values         : "<<outdir<<"/captop_cube_values_f64.raw\n  Occupancy           : "<<outdir<<"/captop_occupied_u8.raw\n  Materials           : "<<outdir<<"/captop_material_i64.raw\n  Element ids         : "<<outdir<<"/captop_element_id_i64.raw\n  Report              : "<<outdir<<"/captop_conversion_report.txt\n\nReady for Stage 5 topology computation: yes\n============================================================\n"; return r.str();
}

static int run_convert(int argc,char** argv){
    if(argc<3){ print_convert_usage(std::cerr); return 1; }
    if(std::string(argv[2])=="--help" || std::string(argv[2])=="-h") { print_convert_usage(std::cout); return 0; }
    std::string input=argv[2]; ParseOptions popts; ValidateOptions vopts; ConvertOptions copts;
    for(int i=3;i<argc;++i){ std::string a=argv[i];
        auto need=[&](const std::string& n){ if(i+1>=argc) throw std::runtime_error(n+" requires a value"); return std::string(argv[++i]); };
        try{
        if(a=="--grid"){ auto v=to_lower(need(a)); if(v=="strict"||v=="strict-cube") vopts.grid_mode=GridMode::StrictCube; else if(v=="rectilinear") vopts.grid_mode=GridMode::Rectilinear; else if(v=="approximate-cube") vopts.grid_mode=GridMode::ApproximateCube; else throw std::runtime_error("unsupported grid mode '"+v+"'"); }
        else if(a=="--tol"){ if(!parse_double(need(a),vopts.tol)||vopts.tol<=0) throw std::runtime_error("invalid --tol"); }
        else if(a=="--cube-rel-tol"){ if(!parse_double(need(a),vopts.cube_rel_tol)||vopts.cube_rel_tol<0) throw std::runtime_error("invalid --cube-rel-tol"); }
        else if(a=="--ignore-non-hexa") popts.ignore_non_hexa=true;
        else if(a=="--max-errors"){ long long n; if(!parse_long_long(need(a),n)||n<=0) throw std::runtime_error("invalid --max-errors"); vopts.max_errors=n; }
        else if(a=="--out") copts.out_dir=need(a);
        else if(a=="--filtration") copts.filtration=need(a);
        else if(a=="--scalar-file") copts.scalar_file=need(a);
        else if(a=="--threshold"){ if(!parse_double(need(a),copts.threshold)) throw std::runtime_error("invalid --threshold"); copts.has_threshold=true; }
        else if(a=="--threshold-op") copts.threshold_op=need(a);
        else if(a=="--selected-material"){ long long m; if(!parse_long_long(need(a),m)) throw std::runtime_error("invalid --selected-material"); copts.has_selected_material=true; copts.selected_material=m; }
        else if(a=="--missing-value"){ if(need(a)!="inf") throw std::runtime_error("--missing-value only accepts inf in Stage 4"); }
        else if(a=="--force") copts.force=true; else if(a=="--overwrite") copts.overwrite=true;
        else if(a=="--max-memory-gb"){ if(!parse_double(need(a),copts.max_memory_gb)||copts.max_memory_gb<=0) throw std::runtime_error("invalid --max-memory-gb"); }
        else throw std::runtime_error("unknown option '"+a+"'");
        } catch(const std::runtime_error& e){ std::cerr<<"error: "<<e.what()<<"\n"; return 1; }
    }
    try{
        if(copts.has_selected_material && copts.filtration!="material") throw std::runtime_error("--selected-material is only compatible with --filtration material");
        if((copts.filtration=="scalar-file"||copts.filtration=="binary-threshold") && copts.scalar_file.empty()) throw std::runtime_error("--filtration "+copts.filtration+" requires --scalar-file");
        if(copts.filtration=="binary-threshold" && !copts.has_threshold) throw std::runtime_error("--filtration binary-threshold requires --threshold");
        if(copts.filtration!="occupancy"&&copts.filtration!="material"&&copts.filtration!="scalar-file"&&copts.filtration!="binary-threshold") throw std::runtime_error("unknown filtration policy '"+copts.filtration+"'");
        Mesh mesh=parse_mesh_file(input,popts); ValidationResult vr=validate_mesh(mesh,vopts); if(!vr.valid){ print_validation_report(input,vopts,vr); return 2; }
        if(vr.nx<=0||vr.ny<=0||vr.nz<=0) throw std::runtime_error("grid dimensions must be positive");
        size_t sx=vr.nx, sy=vr.ny, sz=vr.nz; if(sx>std::numeric_limits<size_t>::max()/sy || sx*sy>std::numeric_limits<size_t>::max()/sz) throw std::runtime_error("integer overflow in nx * ny * nz");
        size_t total=sx*sy*sz, bytes=total*sizeof(GridCell); double limit=copts.max_memory_gb*1024.0*1024.0*1024.0; if(bytes>limit&&!copts.force) throw std::runtime_error("dense memory estimate exceeds --max-memory-gb; use --force to override");
        std::vector<GridCell> cells(total); std::unordered_map<long long,double> scalars; if(copts.filtration=="scalar-file"||copts.filtration=="binary-threshold") scalars=read_scalar_file(copts.scalar_file);
        for(const auto& ic: vr.indexed_cells){ size_t idx=linear_index(ic.i,ic.j,ic.k,vr.nx,vr.ny); auto& g=cells.at(idx); if(g.occupied) throw std::runtime_error("duplicate occupied grid cell "+cell_key_string(CellIndex{ic.i,ic.j,ic.k})); g.occupied=true; g.original_element_id=ic.original_element_id; g.has_material=ic.has_material; g.material=ic.has_material?ic.material:0; double val=0.0;
            if(copts.filtration=="occupancy") val=0.0; else if(copts.filtration=="material"){ if(!ic.has_material) throw std::runtime_error("element "+std::to_string(ic.original_element_id)+" has no material/layer id, but --filtration material requires material values for all occupied cells"); val=copts.has_selected_material ? (ic.material==copts.selected_material?0.0:std::numeric_limits<double>::infinity()) : static_cast<double>(ic.material); }
            else { auto it=scalars.find(ic.original_element_id); if(it==scalars.end()) throw std::runtime_error("missing scalar value for element "+std::to_string(ic.original_element_id)); val=(copts.filtration=="scalar-file")?it->second:(threshold_pass(it->second,copts.threshold,copts.threshold_op,vopts.tol)?0.0:std::numeric_limits<double>::infinity()); }
            g.cube_value=val; }
        ConvertCounts cnt; std::vector<double> values(total); std::vector<uint8_t> occ(total); std::vector<int64_t> mat(total), eid(total);
        for(size_t i=0;i<total;++i){ const auto& g=cells[i]; values[i]=g.cube_value; occ[i]=g.occupied?1:0; mat[i]=g.occupied?(g.has_material?g.material:0):-1; eid[i]=g.occupied?g.original_element_id:-1; if(g.occupied) cnt.occupied++; else cnt.missing++; if(std::isfinite(g.cube_value)) cnt.finite++; }
        cnt.selected=cnt.finite; if(cnt.occupied!=mesh.hexes.size()) throw std::runtime_error("occupied_cubes does not match parsed hexahedra");
        std::filesystem::path outdir(copts.out_dir); std::filesystem::create_directories(outdir); if(!std::filesystem::is_directory(outdir)) throw std::runtime_error("output directory cannot be created");
        std::vector<std::string> names={"captop_grid_metadata.json","captop_cube_values_f64.raw","captop_occupied_u8.raw","captop_material_i64.raw","captop_element_id_i64.raw","captop_conversion_report.txt"}; for(auto& n:names) if(!copts.overwrite && std::filesystem::exists(outdir/n)) throw std::runtime_error("output file already exists: "+(outdir/n).string());
        write_raw(outdir/"captop_cube_values_f64.raw",values); write_raw(outdir/"captop_occupied_u8.raw",occ); write_raw(outdir/"captop_material_i64.raw",mat); write_raw(outdir/"captop_element_id_i64.raw",eid);
        { std::ofstream js(outdir/"captop_grid_metadata.json"); js<<std::setprecision(17)<<"{\n  \"software\": {\"name\": \"captop\", \"version\": \""<<CAPTOP_VERSION<<"\"},\n  \"input\": {\"path\": \""<<json_escape(input)<<"\"},\n  \"validation\": {\"grid_mode\": \""<<grid_mode_name(vopts.grid_mode)<<"\", \"tolerance\": "<<vopts.tol<<", \"status\": \"VALID\"},\n  \"grid\": {\"nx\": "<<vr.nx<<", \"ny\": "<<vr.ny<<", \"nz\": "<<vr.nz<<", \"total_cells\": "<<total<<", \"origin\": ["<<vr.origin.x<<", "<<vr.origin.y<<", "<<vr.origin.z<<"], \"spacing\": ["<<vr.spacing.x<<", "<<vr.spacing.y<<", "<<vr.spacing.z<<"], \"index_order\": \"i + nx * (j + ny * k)\", \"axis_order\": [\"x\", \"y\", \"z\"]}"; if(vopts.grid_mode==GridMode::Rectilinear){ js<<",\n  \"axes\": {\"x\": "; write_axis_json(js,vr.x_axis); js<<", \"y\": "; write_axis_json(js,vr.y_axis); js<<", \"z\": "; write_axis_json(js,vr.z_axis); js<<"}";} js<<",\n  \"counts\": {\"nodes\": "<<vr.n_nodes<<", \"hexahedra\": "<<vr.n_hexes<<", \"occupied_cubes\": "<<cnt.occupied<<", \"missing_cubes\": "<<cnt.missing<<", \"finite_value_cubes\": "<<cnt.finite<<", \"selected_cubes\": "<<cnt.selected<<"},\n  \"filtration\": {\"policy\": \""<<copts.filtration<<"\", \"missing_value\": \"+inf\", \"scalar_file\": "; if(copts.scalar_file.empty()) js<<"null"; else js<<"\""<<json_escape(copts.scalar_file)<<"\""; js<<", \"threshold\": "; if(copts.has_threshold) js<<copts.threshold; else js<<"null"; js<<", \"threshold_op\": "; if(copts.filtration=="binary-threshold") js<<"\""<<copts.threshold_op<<"\""; else js<<"null"; js<<", \"selected_material\": "; if(copts.has_selected_material) js<<copts.selected_material; else js<<"null"; js<<"},\n  \"raw_files\": {\"cube_values_f64\": {\"path\": \"captop_cube_values_f64.raw\", \"type\": \"float64\", \"endianness\": \"little\", \"count\": "<<total<<"}, \"occupied_u8\": {\"path\": \"captop_occupied_u8.raw\", \"type\": \"uint8\", \"count\": "<<total<<"}, \"material_i64\": {\"path\": \"captop_material_i64.raw\", \"type\": \"int64\", \"endianness\": \"little\", \"count\": "<<total<<"}, \"element_id_i64\": {\"path\": \"captop_element_id_i64.raw\", \"type\": \"int64\", \"endianness\": \"little\", \"count\": "<<total<<"}}\n}\n"; if(!js) throw std::runtime_error("cannot write metadata"); }
        std::string rep=conversion_report(input,copts.out_dir,vopts,copts,vr,cnt,bytes); { std::ofstream rr(outdir/"captop_conversion_report.txt"); rr<<rep; if(!rr) throw std::runtime_error("cannot write conversion report"); } std::cout<<rep; return 0;
    } catch(const ParseError& e){ std::cerr<<"parse error: "<<e.what()<<"\n"; return 1; } catch(const std::exception& e){ std::cerr<<"conversion error: "<<e.what()<<"\n"; return 3; }
}

enum class TopologySource { NONE, GUDHI, FALLBACK_DIAGNOSTIC };
static std::string topology_source_value(TopologySource s){ if(s==TopologySource::GUDHI) return "gudhi"; if(s==TopologySource::FALLBACK_DIAGNOSTIC) return "fallback_diagnostic"; return "none"; }
static std::string json_string_array(const std::vector<std::string>& xs){ std::ostringstream o; o<<"["; for(size_t i=0;i<xs.size();++i){ if(i) o<<", "; o<<"\""<<json_escape(xs[i])<<"\""; } o<<"]"; return o.str(); }
static std::string csv_escape(std::string v){ bool q=v.find_first_of(",\n\r\"")!=std::string::npos; size_t pos=0; while((pos=v.find('"',pos))!=std::string::npos){ v.insert(pos,"\""); pos+=2; } return q?"\""+v+"\"":v; }


static bool is_prime_field(int p) {
  if (p < 2) return false;
  for (int d = 2; d * d <= p; ++d) if (p % d == 0) return false;
  return true;
}
struct PersistOptions {
  std::string out_dir="captop_persist_out", filtration="occupancy", scalar_file, mode="sublevel", betti_curve_values="unique";
  std::vector<int> dims{0,1,2};
  int field=2, betti_curve_samples=200; double min_persistence=0, max_memory_gb=4.0;
  bool overwrite=false, force=false, quiet=false, write_pairs=false, write_diagrams=false, write_betti_curve=false, write_barcode_summary=false, write_json=false, write_all=false, finite_only=false, include_essential=true;
};
struct PersistenceInterval { int dim=0; double bc=0, dc=0, bp=0, dp=0, pers=0; bool inf=false; };
static std::string fnum(double v){ if(std::isinf(v)) return v>0?"inf":"-inf"; std::ostringstream o; o<<std::setprecision(17)<<v; return o.str(); }
static std::vector<int> parse_dims(const std::string& s){ std::vector<int> r; std::stringstream ss(s); std::string t; while(std::getline(ss,t,',')){ long long d; if(!parse_long_long(trim(t),d)) throw std::runtime_error("invalid homology dimension '"+t+"'"); if(d<0) throw std::runtime_error("homology dimensions must be nonnegative"); if(d>3) throw std::runtime_error("unsupported homology dimension "+std::to_string(d)+" for 3D cubical grid"); r.push_back((int)d);} if(r.empty()) throw std::runtime_error("--homology-dim cannot be empty"); std::sort(r.begin(),r.end()); r.erase(std::unique(r.begin(),r.end()),r.end()); return r; }
static bool dim_requested(const PersistOptions& o,int d){ return std::find(o.dims.begin(),o.dims.end(),d)!=o.dims.end(); }
static std::string dims_csv(const std::vector<int>& d){ std::ostringstream o; for(size_t i=0;i<d.size();++i){ if(i)o<<","; o<<d[i]; } return o.str(); }
static void print_persist_usage(std::ostream& os){ os<<"Usage:\n  captop persist <input-mesh> [options]\n\nStage 6 persistent homology using GUDHI bitmap cubical complexes.\n\n"; print_validation_options(os); os<<"\n"; print_persist_options(os); os<<"  -h, --help               Show this help message\n"; }

static int run_persist(int argc,char** argv){
  if(argc<3){ print_persist_usage(std::cerr); return 1; }
  if(std::string(argv[2])=="--help"||std::string(argv[2])=="-h"){ print_persist_usage(std::cout); return 0; }
#ifndef CAPTOP_WITH_GUDHI
  std::cerr<<"CAPTOP persistent homology requires GUDHI support.\nRebuild with -DCAPTOP_WITH_GUDHI=ON.\n"; return 1;
#endif
  std::string input=argv[2]; ParseOptions po; ValidateOptions vo; PersistOptions opt;
  for(int i=3;i<argc;++i){ std::string a=argv[i]; auto need=[&](const std::string& n){ if(i+1>=argc) throw std::runtime_error(n+" requires a value"); return std::string(argv[++i]);};
    try{ if(a=="--grid"){ auto v=to_lower(need(a)); if(v=="strict"||v=="strict-cube") vo.grid_mode=GridMode::StrictCube; else if(v=="rectilinear") vo.grid_mode=GridMode::Rectilinear; else throw std::runtime_error("unsupported grid mode"); }
      else if(a=="--tol"){ if(!parse_double(need(a),vo.tol)||vo.tol<=0) throw std::runtime_error("invalid --tol"); }
      else if(a=="--ignore-non-hexa") po.ignore_non_hexa=true; else if(a=="--max-errors"){ long long n; if(!parse_long_long(need(a),n)||n<=0) throw std::runtime_error("invalid --max-errors"); vo.max_errors=n; }
      else if(a=="--homology-dim") opt.dims=parse_dims(need(a)); else if(a=="--field"){ long long f; if(!parse_long_long(need(a),f)||f>INT32_MAX) throw std::runtime_error("invalid --field"); opt.field=(int)f; }
      else if(a=="--filtration"){ opt.filtration=need(a); if(opt.filtration=="scalar-file" && i+1<argc && std::string(argv[i+1]).rfind("--",0)!=0) opt.scalar_file=argv[++i]; }
      else if(a=="--scalar-file") opt.scalar_file=need(a); else if(a=="--mode") opt.mode=to_lower(need(a)); else if(a=="--min-persistence"){ if(!parse_double(need(a),opt.min_persistence)||opt.min_persistence<0) throw std::runtime_error("invalid minimum persistence"); }
      else if(a=="--out") opt.out_dir=need(a); else if(a=="--overwrite") opt.overwrite=true; else if(a=="--force") opt.force=true; else if(a=="--quiet") opt.quiet=true; else if(a=="--max-memory-gb"){ if(!parse_double(need(a),opt.max_memory_gb)||opt.max_memory_gb<=0) throw std::runtime_error("invalid --max-memory-gb"); }
      else if(a=="--write-pairs") opt.write_pairs=true; else if(a=="--write-diagrams") opt.write_diagrams=true; else if(a=="--write-betti-curve") opt.write_betti_curve=true; else if(a=="--write-barcode-summary") opt.write_barcode_summary=true; else if(a=="--write-json") opt.write_json=true; else if(a=="--write-all") opt.write_all=true; else if(a=="--finite-only"){ opt.finite_only=true; opt.include_essential=false; } else if(a=="--include-essential") opt.include_essential=true; else if(a=="--essential-death") (void)need(a);
      else if(a=="--betti-curve-values") opt.betti_curve_values=need(a); else if(a=="--betti-curve-samples"){ long long n; if(!parse_long_long(need(a),n)||n<2) throw std::runtime_error("--betti-curve-samples must be >= 2"); opt.betti_curve_samples=(int)n; } else throw std::runtime_error("unknown option '"+a+"'");
    }catch(const std::exception& e){ std::cerr<<"error: "<<e.what()<<"\n"; return 1; }}
  if(!is_prime_field(opt.field)){ std::cerr<<"error: coefficient field must be prime (got "<<opt.field<<")\n"; return 1; }
  if(opt.mode!="sublevel"&&opt.mode!="superlevel"){ std::cerr<<"error: --mode must be sublevel or superlevel\n"; return 1; }
  if(opt.filtration=="scalar-file"&&opt.scalar_file.empty()){ std::cerr<<"error: --filtration scalar-file requires --scalar-file or inline CSV path\n"; return 3; }
  if(opt.filtration!="occupancy"&&opt.filtration!="material"&&opt.filtration!="scalar-file"){ std::cerr<<"error: unsupported filtration policy '"<<opt.filtration<<"'\n"; return 1; }
  try{
    auto t0=std::chrono::steady_clock::now(); Mesh mesh=parse_mesh_file(input,po); ValidationResult vr=validate_mesh(mesh,vo); if(!vr.valid){ print_validation_report(input,vo,vr); return 2; }
    size_t total=(size_t)vr.nx*(size_t)vr.ny*(size_t)vr.nz; if(total*sizeof(double)>opt.max_memory_gb*1024.0*1024*1024&&!opt.force) throw std::runtime_error("dense grid memory estimate exceeds --max-memory-gb; use --force to override");
    std::vector<double> phys(total,std::numeric_limits<double>::infinity()), comp(total,std::numeric_limits<double>::infinity()); std::unordered_map<long long,double> scalars; if(opt.filtration=="scalar-file") scalars=read_scalar_file(opt.scalar_file);
    for(const auto& ic: vr.indexed_cells){ size_t idx=linear_index(ic.i,ic.j,ic.k,vr.nx,vr.ny); double v=0; if(opt.filtration=="occupancy") v=0; else if(opt.filtration=="material"){ if(!ic.has_material) throw std::runtime_error("element "+std::to_string(ic.original_element_id)+" lacks material value"); v=(double)ic.material; } else { auto it=scalars.find(ic.original_element_id); if(it==scalars.end()) throw std::runtime_error("missing scalar value for element "+std::to_string(ic.original_element_id)); v=it->second; if(!std::isfinite(v)) throw std::runtime_error("scalar value for element "+std::to_string(ic.original_element_id)+" is nonfinite"); } phys[idx]=v; comp[idx]=(opt.mode=="superlevel")?-v:v; }
    std::vector<double> finite; for(double v:phys) if(std::isfinite(v)) finite.push_back(v); double pmin=finite.empty()?0:*std::min_element(finite.begin(),finite.end()), pmax=finite.empty()?0:*std::max_element(finite.begin(),finite.end());
#ifdef CAPTOP_WITH_GUDHI
    using Base=Gudhi::cubical_complex::Bitmap_cubical_complex_base<double>; using Complex=Gudhi::cubical_complex::Bitmap_cubical_complex<Base>; using Field=Gudhi::persistent_cohomology::Field_Zp; using Pcoh=Gudhi::persistent_cohomology::Persistent_cohomology<Complex,Field>;
    std::vector<unsigned> gdims={static_cast<unsigned>(vr.nx),static_cast<unsigned>(vr.ny),static_cast<unsigned>(vr.nz)}; Complex cc(gdims,comp,true); Pcoh pcoh(cc); pcoh.init_coefficients(opt.field); pcoh.compute_persistent_cohomology(opt.min_persistence);
    std::stringstream diag; pcoh.output_diagram(diag); std::vector<PersistenceInterval> ints; std::string line; while(std::getline(diag,line)){ std::istringstream ls(line); int d; double b,death; if(!(ls>>d>>b)) continue; bool inf=false; if(!(ls>>death)){ death=std::numeric_limits<double>::infinity(); inf=true; } if(std::isinf(death)) inf=true; if(!dim_requested(opt,d)) continue; if(opt.finite_only&&inf) continue; double pers=inf?std::numeric_limits<double>::infinity():death-b; if(pers+1e-14<opt.min_persistence) continue; PersistenceInterval pi; pi.dim=d; pi.bc=b; pi.dc=death; pi.inf=inf; pi.pers=pers; pi.bp=(opt.mode=="superlevel"?-b:b); pi.dp=inf?(opt.mode=="superlevel"?-std::numeric_limits<double>::infinity():std::numeric_limits<double>::infinity()):(opt.mode=="superlevel"?-death:death); ints.push_back(pi); }
    std::sort(ints.begin(),ints.end(),[](const auto&a,const auto&b){ return std::make_tuple(a.dim,a.bc,a.inf,a.dc)<std::make_tuple(b.dim,b.bc,b.inf,b.dc); });
    std::filesystem::path od(opt.out_dir); std::filesystem::create_directories(od); std::vector<std::string> names={"persistence_pairs.csv","barcode_summary.csv","betti_curve.csv","captop_persistence_summary.json","captop_persistence_report.txt"}; for(int d:opt.dims) names.push_back("diagram_dim"+std::to_string(d)+".csv"); for(auto&n:names) if(!opt.overwrite&&std::filesystem::exists(od/n)) throw std::runtime_error("output file already exists: "+(od/n).string());
    { std::ofstream f(od/"persistence_pairs.csv"); f<<"pair_id,dimension,birth_computational,death_computational,birth_physical,death_physical,death_type,persistence_computational,persistence_physical_abs,filtration_policy,mode,coefficient_field\n"; int id=0; for(auto&x:ints) f<<id++<<","<<x.dim<<","<<fnum(x.bc)<<","<<fnum(x.dc)<<","<<fnum(x.bp)<<","<<fnum(x.dp)<<","<<(x.inf?"essential":"finite")<<","<<fnum(x.pers)<<","<<fnum(x.inf?std::numeric_limits<double>::infinity():std::abs(x.dp-x.bp))<<","<<opt.filtration<<","<<opt.mode<<","<<opt.field<<"\n"; }
    for(int d:opt.dims){ std::ofstream f(od/("diagram_dim"+std::to_string(d)+".csv")); f<<"birth,death,birth_physical,death_physical,death_type,persistence\n"; for(auto&x:ints) if(x.dim==d) f<<fnum(x.bc)<<","<<fnum(x.dc)<<","<<fnum(x.bp)<<","<<fnum(x.dp)<<","<<(x.inf?"essential":"finite")<<","<<fnum(x.pers)<<"\n"; }
    { std::ofstream f(od/"barcode_summary.csv"); f<<"dimension,intervals_total,intervals_finite,intervals_essential,persistence_min,persistence_q1,persistence_median,persistence_mean,persistence_q3,persistence_max,birth_min,birth_max,death_min,death_max\n"; for(int d:opt.dims){ std::vector<double> ps,bs,ds; long long ess=0; for(auto&x:ints) if(x.dim==d){bs.push_back(x.bc); if(x.inf) ess++; else {ps.push_back(x.pers); ds.push_back(x.dc);}} auto stat=[&](std::vector<double> v,int q){ if(v.empty()) return std::string("NA"); std::sort(v.begin(),v.end()); if(q==0) return fnum(v.front()); if(q==4) return fnum(v.back()); if(q==2) return fnum(v[v.size()/2]); if(q==5) return fnum(std::accumulate(v.begin(),v.end(),0.0)/v.size()); return fnum(v[(v.size()*q)/4]);}; f<<d<<","<<(bs.size())<<","<<ps.size()<<","<<ess<<","<<stat(ps,0)<<","<<stat(ps,1)<<","<<stat(ps,2)<<","<<stat(ps,5)<<","<<stat(ps,3)<<","<<stat(ps,4)<<","<<stat(bs,0)<<","<<stat(bs,4)<<","<<stat(ds,0)<<","<<stat(ds,4)<<"\n"; }}
    std::vector<double> th; for(double v:comp) if(std::isfinite(v)) th.push_back(v); std::sort(th.begin(),th.end()); th.erase(std::unique(th.begin(),th.end()),th.end()); if(opt.betti_curve_values=="uniform"&&!th.empty()){ double a=th.front(),b=th.back(); th.clear(); for(int i=0;i<opt.betti_curve_samples;i++) th.push_back(a+(b-a)*i/(opt.betti_curve_samples-1)); }
    { std::ofstream f(od/"betti_curve.csv"); f<<"threshold_computational,threshold_physical,dimension,betti,mode,filtration_policy\n"; for(double t:th) for(int d:opt.dims){ long long beta=0; for(auto&x:ints) if(x.dim==d && x.bc<=t && (x.inf || t<x.dc)) beta++; f<<fnum(t)<<","<<fnum(opt.mode=="superlevel"?-t:t)<<","<<d<<","<<beta<<","<<opt.mode<<","<<opt.filtration<<"\n"; }}
    long long finite_n=0, ess_n=0; std::map<int,long long> bydim; for(auto&x:ints){ bydim[x.dim]++; if(x.inf) ess_n++; else finite_n++; }
    { std::ofstream j(od/"captop_persistence_summary.json"); j<<std::setprecision(17)<<"{\n  \"software\": {\"name\": \"captop\", \"version\": \""<<CAPTOP_VERSION<<"\"},\n  \"input\": {\"path\": \""<<json_escape(input)<<"\"},\n  \"validation\": {\"grid_mode\": \""<<grid_mode_name(vo.grid_mode)<<"\", \"tolerance\": "<<vo.tol<<", \"status\": \"VALID\"},\n  \"grid\": {\"nx\": "<<vr.nx<<", \"ny\": "<<vr.ny<<", \"nz\": "<<vr.nz<<", \"total_voxels\": "<<total<<", \"occupied_voxels\": "<<vr.indexed_cells.size()<<", \"missing_voxels\": "<<(total-vr.indexed_cells.size())<<", \"occupied_fraction\": "<<(total?double(vr.indexed_cells.size())/total:0)<<", \"origin\": ["<<vr.origin.x<<","<<vr.origin.y<<","<<vr.origin.z<<"], \"spacing\": ["<<vr.spacing.x<<","<<vr.spacing.y<<","<<vr.spacing.z<<"], \"index_order\": \"i + nx * (j + ny * k)\"},\n  \"gudhi\": {\"compiled\": true, \"used\": true, \"success\": true, \"coefficient_field\": "<<opt.field<<", \"input_top_cells\": "<<total<<", \"missing_value\": \"+inf\"},\n  \"filtration\": {\"policy\": \""<<opt.filtration<<"\", \"mode\": \""<<opt.mode<<"\", \"min_persistence\": "<<opt.min_persistence<<", \"scalar_file\": "; if(opt.scalar_file.empty()) j<<"null"; else j<<"\""<<json_escape(opt.scalar_file)<<"\""; j<<", \"finite_value_count\": "<<finite.size()<<", \"infinite_value_count\": "<<(total-finite.size())<<", \"physical_min\": "<<pmin<<", \"physical_max\": "<<pmax<<", \"computational_min\": "<<(th.empty()?0:th.front())<<", \"computational_max\": "<<(th.empty()?0:th.back())<<"},\n  \"persistence\": {\"requested_dimensions\": ["<<dims_csv(opt.dims)<<"], \"interval_count_total\": "<<ints.size()<<", \"interval_count_finite\": "<<finite_n<<", \"interval_count_essential\": "<<ess_n<<", \"intervals_by_dimension\": {"; bool first=true; for(auto&kv:bydim){ if(!first) j<<","; first=false; j<<"\""<<kv.first<<"\":"<<kv.second;} j<<"}},\n  \"outputs\": {\"persistence_pairs_csv\": \""<<(od/"persistence_pairs.csv").string()<<"\", \"barcode_summary_csv\": \""<<(od/"barcode_summary.csv").string()<<"\", \"betti_curve_csv\": \""<<(od/"betti_curve.csv").string()<<"\"},\n  \"timings\": {\"total_seconds\": "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count()<<"},\n  \"status\": \"SUCCESS\"\n}\n"; }
    std::ostringstream rep; rep<<"============================================================\nCAPTOP persistent homology report\n============================================================\nSoftware version      : "<<CAPTOP_VERSION<<"\nInput file            : "<<input<<"\nGrid mode             : "<<grid_mode_name(vo.grid_mode)<<"\nTolerance             : "<<vo.tol<<"\nCoefficient field     : Z/"<<opt.field<<"Z\nFiltration policy     : "<<opt.filtration<<"\nFiltration mode       : "<<opt.mode<<"\nMin persistence       : "<<opt.min_persistence<<"\nPrimary topology      : closed occupied cubical complex\nStatus                : VALID AND GUDHI-PERSISTENCE-ANALYZED\n\nTopology engine\n  GUDHI support       : enabled\n  GUDHI used          : yes\n  Result authority    : GUDHI persistent cohomology\n\nGrid\n  Dimensions          : "<<vr.nx<<" x "<<vr.ny<<" x "<<vr.nz<<"\n  Total voxels        : "<<total<<"\n  Occupied voxels     : "<<vr.indexed_cells.size()<<"\n  Missing voxels      : "<<(total-vr.indexed_cells.size())<<"\n  Occupied fraction   : "<<(total?double(vr.indexed_cells.size())/total:0)<<"\n\nFiltration\n  Policy              : "<<opt.filtration<<"\n  Mode                : "<<opt.mode<<"\n"; if(opt.mode=="superlevel") rep<<"  Superlevel handling   : finite scalar values internally negated; reported values restored to original physical scale\n"; rep<<"  Finite cube values  : "<<finite.size()<<"\n  Infinite values     : "<<(total-finite.size())<<"\n  Physical min/max    : "<<pmin<<" / "<<pmax<<"\n\nPersistence intervals\n  Requested dimensions: "<<dims_csv(opt.dims)<<"\n  Total intervals     : "<<ints.size()<<"\n  Finite intervals    : "<<finite_n<<"\n  Essential intervals : "<<ess_n<<"\n"; for(int d:opt.dims) rep<<"  Dim "<<d<<" intervals     : "<<bydim[d]<<"\n"; rep<<"\nOutput files\n  Pairs               : "<<(od/"persistence_pairs.csv").string()<<"\n  Barcode summary     : "<<(od/"barcode_summary.csv").string()<<"\n  Betti curve         : "<<(od/"betti_curve.csv").string()<<"\n  JSON summary        : "<<(od/"captop_persistence_summary.json").string()<<"\n  Report              : "<<(od/"captop_persistence_report.txt").string()<<"\n\nReady for Stage 7 performance engineering: yes\n============================================================\n";
    { std::ofstream r(od/"captop_persistence_report.txt"); r<<rep.str(); } if(!opt.quiet) std::cout<<rep.str(); return 0;
#endif
  }catch(const ParseError& e){ std::cerr<<"parse error: "<<e.what()<<"\n"; return 1; }catch(const ValidationError& e){ std::cerr<<"validation error: "<<e.what()<<"\n"; return 2; }catch(const std::exception& e){ std::cerr<<"persistence error: "<<e.what()<<"\nStatus                : FAILED\n"; return 3; }
}
struct BettiOptions { std::string out_dir="captop_betti_out"; bool overwrite=false, force=false, write_diagram=false, write_json=false, write_csv=false, quiet=false, require_gudhi=false, allow_fallback=true, write_contact_audit=false; double max_memory_gb=2.0; int field=2; };
enum class ContactType { Face, Edge, Vertex };
static std::string contact_type_name(ContactType t){ if(t==ContactType::Face) return "face"; if(t==ContactType::Edge) return "edge"; return "vertex"; }
struct ContactRecord {
  long long component_face_a = 0, component_face_b = 0;
  ContactType contact_type = ContactType::Face;
  long long i1 = 0, j1 = 0, k1 = 0, i2 = 0, j2 = 0, k2 = 0, element_id_1 = 0,
            element_id_2 = 0, linear_index_1 = 0, linear_index_2 = 0;
};
struct ComponentMergerRecord {
  long long merger_id = 0, face_component_a = 0, face_component_b = 0,
            closed_component_before_a = 0, closed_component_before_b = 0,
            closed_component_after = 0;
  ContactType contact_type = ContactType::Face;
  long long i1 = 0, j1 = 0, k1 = 0, i2 = 0, j2 = 0, k2 = 0, element_id_1 = 0,
            element_id_2 = 0, linear_index_1 = 0, linear_index_2 = 0;
};
struct ContactAuditResult {
  bool computed = false;
  long long H0_face = 0, H0_closed = 0, face_contact_count = 0,
            edge_contact_count = 0, vertex_contact_count = 0,
            inter_face_component_contact_count = 0,
            inter_face_component_face_contact_count = 0,
            inter_face_component_edge_contact_count = 0,
            inter_face_component_vertex_contact_count = 0,
            face_component_merger_count = 0, merger_face_record_count = 0,
            merger_edge_record_count = 0, merger_vertex_record_count = 0,
            component_reduction = 0;
  std::string contact_merger_file = "captop_contact_mergers.csv",
              contact_records_file = "captop_contact_contacts.csv";
  bool contact_merger_file_written = false,
       contact_records_file_written = false;
  std::vector<ContactRecord> inter_component_contacts;
  std::vector<ComponentMergerRecord> merger_records;
  bool gudhi_H0_matches_closed = false, gudhi_H0_matches_face = false;
};
struct BettiResult {
  bool success = false, gudhi_compiled = false, gudhi_used = false,
       gudhi_success = false, fallback_used = true,
       topology_authoritative = false, closed_cube_topology_primary = true,
       H0_closed_check = false, H0_face_check = false;
  TopologySource topology_source = TopologySource::NONE;
  std::string result_authority, readiness_stage6;
  int H0 = 0, H1 = 0, H2 = 0, H3 = 0, H0_union_find = 0, H0_face = 0,
      H0_closed = 0, H0_gudhi = -1, H2_complement = -1;
  long long result_chi = 0, chi_cells = 0, occupied_voxels = 0,
            missing_voxels = 0, total_voxels = 0, surface_voxels = 0,
            surface_faces = 0, C0 = 0, C1 = 0, C2 = 0, C3 = 0;
  double occupied_fraction = 0;
  std::string H0_source, H1_source, H2_source, H3_source, chi_source;
  ContactAuditResult contact_audit;
  std::vector<std::string> warnings, errors;
};
static bool is_prime_field_betti(int p) {
  if (p < 2)
    return false;
  for (int d = 2; d * d <= p; ++d)
    if (p % d == 0)
      return false;
  return true;
}
struct DSU {
  std::vector<int> p, r;
  explicit DSU(size_t n) : p(n), r(n, 0) { std::iota(p.begin(), p.end(), 0); }
  int find(int x) { return p[x] == x ? x : p[x] = find(p[x]); }
  bool unite(int a, int b) {
    a = find(a);
    b = find(b);
    if (a == b)
      return false;
    if (r[a] < r[b])
      std::swap(a, b);
    p[b] = a;
    if (r[a] == r[b])
      r[a]++;
    return true;
  }
};
static BettiResult analyze_betti(const ValidationResult &vr,
                                 const BettiOptions &opts) {
  BettiResult br;
#ifdef CAPTOP_WITH_GUDHI
  br.gudhi_compiled = true;
#endif
  if (!br.gudhi_compiled && opts.require_gudhi)
    throw std::runtime_error("CAPTOP was built without GUDHI support; rebuild "
                             "with -DCAPTOP_WITH_GUDHI=ON or omit "
                             "--require-gudhi to allow fallback diagnostics.");
  if (vr.nx <= 0 || vr.ny <= 0 || vr.nz <= 0)
    throw std::runtime_error("grid dimensions must be positive");
  size_t sx = vr.nx, sy = vr.ny, sz = vr.nz;
  if (sx > SIZE_MAX / sy || sx * sy > SIZE_MAX / sz)
    throw std::runtime_error("integer overflow in nx * ny * nz");
  size_t total = sx * sy * sz;
  size_t bytes = total * (sizeof(uint8_t) + sizeof(int));
  if (bytes > opts.max_memory_gb * 1024.0 * 1024 * 1024 && !opts.force)
    throw std::runtime_error("dense memory estimate exceeds --max-memory-gb; "
                             "use --force to override");
  br.total_voxels = total;
  br.occupied_voxels = vr.indexed_cells.size();
  br.missing_voxels = br.total_voxels - br.occupied_voxels;
  br.occupied_fraction = total ? double(br.occupied_voxels) / double(total) : 0;
  std::vector<uint8_t> occ(total, 0);
  std::vector<int> compact(total, -1);
  std::vector<long long> element_id(total, 0);
  int cid = 0;
  for (auto &ic : vr.indexed_cells) {
    size_t idx = linear_index(ic.i, ic.j, ic.k, vr.nx, vr.ny);
    if (occ[idx])
      throw std::runtime_error(
          "duplicate occupied grid cell during Betti grid construction");
    occ[idx] = 1;
    compact[idx] = cid++;
    element_id[idx] = ic.original_element_id;
  }
  auto occupied = [&](long long i, long long j, long long k) {
    return i >= 0 && j >= 0 && k >= 0 && i < vr.nx && j < vr.ny && k < vr.nz &&
           occ[linear_index(i, j, k, vr.nx, vr.ny)];
  };
  DSU face_dsu(br.occupied_voxels), closed_dsu(br.occupied_voxels);
  for (auto &ic : vr.indexed_cells) {
    size_t a = linear_index(ic.i, ic.j, ic.k, vr.nx, vr.ny);
    bool surf = false;
    const int dirs[6][3] = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                            {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};
    for (auto &d : dirs) {
      long long ni = ic.i + d[0], nj = ic.j + d[1], nk = ic.k + d[2];
      if (!occupied(ni, nj, nk)) {
        br.surface_faces++;
        surf = true;
      }
    }
    if (surf)
      br.surface_voxels++;
    const int pd[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (auto &d : pd) {
      long long ni = ic.i + d[0], nj = ic.j + d[1], nk = ic.k + d[2];
      if (occupied(ni, nj, nk))
        face_dsu.unite(compact[a],
                       compact[linear_index(ni, nj, nk, vr.nx, vr.ny)]);
    }
  }
  std::map<int, long long> root_min;
  for (auto &ic : vr.indexed_cells) {
    size_t idx = linear_index(ic.i, ic.j, ic.k, vr.nx, vr.ny);
    int r = face_dsu.find(compact[idx]);
    if (!root_min.count(r) || static_cast<long long>(idx) < root_min[r])
      root_min[r] = idx;
  }
  std::vector<std::pair<long long, int>> mins;
  for (auto &kv : root_min)
    mins.push_back({kv.second, kv.first});
  std::sort(mins.begin(), mins.end());
  std::unordered_map<int, int> root_to_face;
  for (size_t i = 0; i < mins.size(); ++i)
    root_to_face[mins[i].second] = static_cast<int>(i);
  std::vector<int> face_comp(total, -1);
  for (auto &ic : vr.indexed_cells) {
    size_t idx = linear_index(ic.i, ic.j, ic.k, vr.nx, vr.ny);
    face_comp[idx] = root_to_face[face_dsu.find(compact[idx])];
  }
  br.H0_face = br.H0_union_find = static_cast<int>(mins.size());
  ContactAuditResult audit;
  audit.computed = true;
  audit.H0_face = br.H0_face;
  for (auto &ic : vr.indexed_cells) {
    size_t aidx = linear_index(ic.i, ic.j, ic.k, vr.nx, vr.ny);
    for (int di = -1; di <= 1; ++di)
      for (int dj = -1; dj <= 1; ++dj)
        for (int dk = -1; dk <= 1; ++dk) {
          if (di == 0 && dj == 0 && dk == 0)
            continue;
          long long ni = ic.i + di, nj = ic.j + dj, nk = ic.k + dk;
          if (!occupied(ni, nj, nk))
            continue;
          size_t bidx = linear_index(ni, nj, nk, vr.nx, vr.ny);
          if (bidx <= aidx)
            continue;
          int m = std::abs(di) + std::abs(dj) + std::abs(dk);
          ContactType ct =
              m == 1 ? ContactType::Face
                     : (m == 2 ? ContactType::Edge : ContactType::Vertex);
          if (m < 1 || m > 3)
            throw std::runtime_error("invalid contact classification");
          if (ct == ContactType::Face)
            audit.face_contact_count++;
          else if (ct == ContactType::Edge)
            audit.edge_contact_count++;
          else
            audit.vertex_contact_count++;
          closed_dsu.unite(compact[aidx], compact[bidx]);
          int fa = face_comp[aidx], fb = face_comp[bidx];
          if (fa != fb) {
            if (fa > fb)
              std::swap(fa, fb);
            ContactRecord cr;
            cr.component_face_a = fa;
            cr.component_face_b = fb;
            cr.contact_type = ct;
            cr.i1 = ic.i;
            cr.j1 = ic.j;
            cr.k1 = ic.k;
            cr.i2 = ni;
            cr.j2 = nj;
            cr.k2 = nk;
            cr.element_id_1 = element_id[aidx];
            cr.element_id_2 = element_id[bidx];
            cr.linear_index_1 = aidx;
            cr.linear_index_2 = bidx;
            audit.inter_component_contacts.push_back(cr);
            audit.inter_face_component_contact_count++;
            if (ct == ContactType::Face) {
              audit.inter_face_component_face_contact_count++;
              br.errors.push_back(
                  "face contact found across distinct face components");
            } else if (ct == ContactType::Edge)
              audit.inter_face_component_edge_contact_count++;
            else
              audit.inter_face_component_vertex_contact_count++;
          }
        }
  }
  std::unordered_set<int> closed_roots;
  for (auto &ic : vr.indexed_cells) {
    size_t idx = linear_index(ic.i, ic.j, ic.k, vr.nx, vr.ny);
    closed_roots.insert(closed_dsu.find(compact[idx]));
  }
  br.H0_closed = static_cast<int>(closed_roots.size());
  audit.H0_closed = br.H0_closed;
  audit.component_reduction = br.H0_face - br.H0_closed;
  auto pri = [](ContactType t) {
    return t == ContactType::Edge ? 0 : (t == ContactType::Vertex ? 1 : 2);
  };
  std::sort(audit.inter_component_contacts.begin(),
            audit.inter_component_contacts.end(),
            [&](const ContactRecord &a, const ContactRecord &b) {
              return std::make_tuple(pri(a.contact_type), a.component_face_a,
                                     a.component_face_b, a.linear_index_1,
                                     a.linear_index_2) <
                     std::make_tuple(pri(b.contact_type), b.component_face_a,
                                     b.component_face_b, b.linear_index_1,
                                     b.linear_index_2);
            });
  DSU comp_dsu(br.H0_face);
  long long mid = 1;
  for (const auto &cr : audit.inter_component_contacts) {
    int ra = comp_dsu.find(cr.component_face_a),
        rb = comp_dsu.find(cr.component_face_b);
    if (ra == rb)
      continue;
    ComponentMergerRecord mr;
    mr.merger_id = mid++;
    mr.face_component_a = cr.component_face_a;
    mr.face_component_b = cr.component_face_b;
    mr.closed_component_before_a = ra;
    mr.closed_component_before_b = rb;
    comp_dsu.unite(ra, rb);
    mr.closed_component_after = comp_dsu.find(ra);
    mr.contact_type = cr.contact_type;
    mr.i1 = cr.i1;
    mr.j1 = cr.j1;
    mr.k1 = cr.k1;
    mr.i2 = cr.i2;
    mr.j2 = cr.j2;
    mr.k2 = cr.k2;
    mr.element_id_1 = cr.element_id_1;
    mr.element_id_2 = cr.element_id_2;
    mr.linear_index_1 = cr.linear_index_1;
    mr.linear_index_2 = cr.linear_index_2;
    audit.merger_records.push_back(mr);
  }
  audit.face_component_merger_count = audit.merger_records.size();
  for (const auto &mr : audit.merger_records) {
    if (mr.contact_type == ContactType::Face)
      audit.merger_face_record_count++;
    else if (mr.contact_type == ContactType::Edge)
      audit.merger_edge_record_count++;
    else
      audit.merger_vertex_record_count++;
  }
  br.contact_audit = audit;
  std::unordered_set<std::string> V, E, F;
  for (auto &ic : vr.indexed_cells) {
    long long i = ic.i, j = ic.j, k = ic.k;
    for (int dx = 0; dx < 2; ++dx)
      for (int dy = 0; dy < 2; ++dy)
        for (int dz = 0; dz < 2; ++dz)
          V.insert(std::to_string(i + dx) + "," + std::to_string(j + dy) + "," +
                   std::to_string(k + dz));
    for (int dy = 0; dy < 2; ++dy)
      for (int dz = 0; dz < 2; ++dz)
        E.insert("x," + std::to_string(i) + "," + std::to_string(j + dy) + "," +
                 std::to_string(k + dz));
    for (int dx = 0; dx < 2; ++dx)
      for (int dz = 0; dz < 2; ++dz)
        E.insert("y," + std::to_string(i + dx) + "," + std::to_string(j) + "," +
                 std::to_string(k + dz));
    for (int dx = 0; dx < 2; ++dx)
      for (int dy = 0; dy < 2; ++dy)
        E.insert("z," + std::to_string(i + dx) + "," + std::to_string(j + dy) +
                 "," + std::to_string(k));
    for (int dx = 0; dx < 2; ++dx)
      F.insert("x," + std::to_string(i + dx) + "," + std::to_string(j) + "," +
               std::to_string(k));
    for (int dy = 0; dy < 2; ++dy)
      F.insert("y," + std::to_string(i) + "," + std::to_string(j + dy) + "," +
               std::to_string(k));
    for (int dz = 0; dz < 2; ++dz)
      F.insert("z," + std::to_string(i) + "," + std::to_string(j) + "," +
               std::to_string(k + dz));
  }
  br.C0 = V.size();
  br.C1 = E.size();
  br.C2 = F.size();
  br.C3 = br.occupied_voxels;
  br.chi_cells = br.C0 - br.C1 + br.C2 - br.C3;
  long long px = vr.nx + 2, py = vr.ny + 2, pz = vr.nz + 2;
  std::vector<uint8_t> block(px * py * pz, 0), seen(px * py * pz, 0);
  auto pidx = [&](long long i, long long j, long long k) {
    return size_t(i + px * (j + py * k));
  };
  for (auto &ic : vr.indexed_cells)
    block[pidx(ic.i + 1, ic.j + 1, ic.k + 1)] = 1;
  std::queue<std::array<long long, 3>> q;
  q.push({0, 0, 0});
  seen[pidx(0, 0, 0)] = 1;
  const int dirs[6][3] = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                          {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};
  while (!q.empty()) {
    auto a = q.front();
    q.pop();
    for (auto &d : dirs) {
      long long ni = a[0] + d[0], nj = a[1] + d[1], nk = a[2] + d[2];
      if (ni >= 0 && nj >= 0 && nk >= 0 && ni < px && nj < py && nk < pz &&
          !block[pidx(ni, nj, nk)] && !seen[pidx(ni, nj, nk)]) {
        seen[pidx(ni, nj, nk)] = 1;
        q.push({ni, nj, nk});
      }
    }
  }
  int cavities = 0;
  for (long long k = 1; k <= vr.nz; ++k)
    for (long long j = 1; j <= vr.ny; ++j)
      for (long long i = 1; i <= vr.nx; ++i)
        if (!block[pidx(i, j, k)] && !seen[pidx(i, j, k)]) {
          cavities++;
          seen[pidx(i, j, k)] = 1;
          q.push({i, j, k});
          while (!q.empty()) {
            auto a = q.front();
            q.pop();
            for (auto &d : dirs) {
              long long ni = a[0] + d[0], nj = a[1] + d[1], nk = a[2] + d[2];
              if (ni >= 1 && nj >= 1 && nk >= 1 && ni <= vr.nx && nj <= vr.ny &&
                  nk <= vr.nz && !block[pidx(ni, nj, nk)] &&
                  !seen[pidx(ni, nj, nk)]) {
                seen[pidx(ni, nj, nk)] = 1;
                q.push({ni, nj, nk});
              }
            }
          }
        }
  br.H2_complement = cavities;
#ifdef CAPTOP_WITH_GUDHI
  using Base = Gudhi::cubical_complex::Bitmap_cubical_complex_base<double>;
  using Complex = Gudhi::cubical_complex::Bitmap_cubical_complex<Base>;
  using Field = Gudhi::persistent_cohomology::Field_Zp;
  using Pcoh =
      Gudhi::persistent_cohomology::Persistent_cohomology<Complex, Field>;
  std::vector<unsigned> dims = {static_cast<unsigned>(vr.nx),
                                static_cast<unsigned>(vr.ny),
                                static_cast<unsigned>(vr.nz)};
  std::vector<double> top(total, std::numeric_limits<double>::infinity());
  for (size_t i = 0; i < total; ++i)
    if (occ[i])
      top[i] = 0.0;
  Complex cc(dims, top, true);
  Pcoh pcoh(cc);
  pcoh.init_coefficients(opts.field);
  pcoh.compute_persistent_cohomology(-1.0);
  br.H0 = pcoh.persistent_betti_number(0, 0.0, 0.0);
  br.H0_gudhi = br.H0;
  br.H1 = pcoh.persistent_betti_number(1, 0.0, 0.0);
  br.H2 = pcoh.persistent_betti_number(2, 0.0, 0.0);
  br.H3 = pcoh.persistent_betti_number(3, 0.0, 0.0);
  br.gudhi_used = br.gudhi_success = true;
  br.topology_source = TopologySource::GUDHI;
  br.topology_authoritative = true;
  br.result_authority = "GUDHI persistent Betti numbers";
  br.H0_source = br.H1_source = br.H2_source = br.H3_source =
      "gudhi_persistent_betti";
  br.chi_source = "gudhi_betti_numbers";
#else
  br.warnings.push_back(
      "CAPTOP was built without GUDHI support. Betti numbers are diagnostic "
      "fallback estimates, not GUDHI-authenticated results.");
  br.H0 = br.H0_closed;
  br.H2 = br.H2_complement;
  br.H3 = 0;
  br.H1 = br.H0_closed + br.H2 - br.chi_cells;
  br.topology_source = TopologySource::FALLBACK_DIAGNOSTIC;
  br.topology_authoritative = false;
  br.result_authority = "diagnostic fallback, not GUDHI-authenticated";
  br.H0_source = "closed_cube_26_neighbor_union_find";
  br.H1_source = "inferred_from_closed_H0_complement_H2_and_cellular_euler";
  br.H2_source = "complement_flood_fill";
  br.H3_source = "assumed_zero_fallback";
  br.chi_source = "cubical_cell_count";
#endif
  br.result_chi = br.H0 - br.H1 + br.H2 - br.H3;
  br.H0_closed_check = br.gudhi_used ? (br.H0_gudhi == br.H0_closed) : false;
  br.H0_face_check = br.gudhi_used ? (br.H0_gudhi == br.H0_face) : false;
  br.contact_audit.gudhi_H0_matches_closed = br.H0_closed_check;
  br.contact_audit.gudhi_H0_matches_face = br.H0_face_check;
  if (br.H0_closed > br.H0_face)
    br.errors.push_back("contact audit impossible: H0_closed exceeds H0_face");
  if (br.contact_audit.component_reduction != br.H0_face - br.H0_closed)
    br.errors.push_back("contact audit component reduction is inconsistent");
  if (br.contact_audit.inter_face_component_contact_count !=
      br.contact_audit.inter_face_component_face_contact_count +
          br.contact_audit.inter_face_component_edge_contact_count +
          br.contact_audit.inter_face_component_vertex_contact_count)
    br.errors.push_back(
        "contact audit raw contact decomposition is inconsistent");
  if (static_cast<long long>(br.contact_audit.merger_records.size()) !=
      br.contact_audit.component_reduction)
    br.errors.push_back("contact audit merger record count is inconsistent");
  if (br.contact_audit.face_component_merger_count !=
      br.contact_audit.merger_face_record_count +
          br.contact_audit.merger_edge_record_count +
          br.contact_audit.merger_vertex_record_count)
    br.errors.push_back("contact audit merger decomposition is inconsistent");
  if (br.contact_audit.merger_face_record_count > 0)
    br.errors.push_back("contact audit selected a face-contact merger between "
                        "distinct face components");
  if (br.H0_face > br.H0_closed &&
      br.contact_audit.inter_face_component_edge_contact_count +
              br.contact_audit.inter_face_component_vertex_contact_count <=
          0)
    br.errors.push_back("contact audit does not explain H0 reduction with "
                        "edge/vertex contacts");
  if (br.gudhi_used) {
    if (br.H0_gudhi != br.H0_closed)
      br.errors.push_back("GUDHI H0 disagrees with closed-cube H0");
    if (br.result_chi != br.chi_cells)
      br.errors.push_back("GUDHI Euler characteristic disagrees with "
                          "cubical-cell Euler characteristic");
  } else {
    if (br.result_chi != br.chi_cells)
      br.errors.push_back(
          "fallback Euler inference is internally inconsistent");
  }
  br.readiness_stage6 =
      br.gudhi_compiled ? "yes" : "no, rebuild with -DCAPTOP_WITH_GUDHI=ON";
  br.success = br.errors.empty();
  return br;
}
static std::string h0_face_note(const BettiResult &br) {
  long long d = br.H0_face - br.H0_closed;
  if (d == 0)
    return "same as closed-cube H0";
  std::ostringstream o;
  o << "differs by " << d << " due to edge/vertex contacts";
  return o.str();
}
static bool should_write_contact_audit_files(const BettiOptions &bo,
                                             const BettiResult &br) {
  return bo.write_contact_audit || br.H0_face != br.H0_closed ||
         (br.gudhi_used && br.H0_gudhi != br.H0_closed) || !br.errors.empty();
}
static std::string contact_audit_file_label(bool written,
                                            const std::string &name) {
  return written ? name : "not written";
}
static std::string betti_report(const std::string &input,
                                const ValidateOptions &vo,
                                const BettiOptions &bo,
                                const ValidationResult &vr,
                                const BettiResult &br) {
  std::ostringstream o;
  bool fb = br.topology_source == TopologySource::FALLBACK_DIAGNOSTIC;
  bool write_contact_files = should_write_contact_audit_files(bo, br);
  o << "============================================================\nCAPTOP "
       "Betti topology "
       "report\n============================================================\n"
    << "Software version      : " << CAPTOP_VERSION
    << "\nInput file            : " << input
    << "\nGrid mode             : " << grid_mode_name(vo.grid_mode)
    << "\nTolerance             : " << std::setprecision(12) << vo.tol << "\n";
  if (fb)
    o << "Requested field       : Z/" << bo.field
      << "Z\nCoefficient field     : not used; GUDHI disabled\n";
  else
    o << "Coefficient field     : Z/" << bo.field << "Z\n";
  o << "Primary topology      : closed occupied cubical complex\nPrimary "
       "connectivity  : 26-neighbor closed-cube connectivity\nFace-adj. "
       "diagnostic  : 6-neighbor voxel connectivity\nStatus                : "
    << (br.success ? (fb ? "VALID AND DIAGNOSTICALLY ANALYZED"
                         : "VALID AND GUDHI-ANALYZED")
                   : "FAILED")
    << "\n";
  for (const auto &w : br.warnings)
    o << "Warning               : " << w << "\n";
  o << "\nTopology engine\n  GUDHI support       : "
    << (br.gudhi_compiled ? "enabled" : "disabled at compile time")
    << "\n  GUDHI used          : " << (br.gudhi_used ? "yes" : "no")
    << "\n  Fallback used       : "
    << (br.gudhi_used ? "yes, for cross-checks" : "yes")
    << "\n  Result authority    : " << br.result_authority
    << "\n\nGrid\n  Dimensions          : " << vr.nx << " x " << vr.ny << " x "
    << vr.nz << "\n  Total voxels        : " << br.total_voxels
    << "\n  Occupied voxels     : " << br.occupied_voxels
    << "\n  Missing voxels      : " << br.missing_voxels
    << "\n  Occupied fraction   : " << br.occupied_fraction << "\n\n";
  if (fb)
    o << "Topology inferred from fallback diagnostics\n  H0 closed-cube "
         "components : "
      << br.H0
      << "  [26-neighbor]\n  H0 face components        : " << br.H0_face
      << "  [6-neighbor diagnostic]\n  H1 tunnels                : " << br.H1
      << "  [inferred from H0_closed + H2 - chi]\n  H2 cavities               "
         ": "
      << br.H2
      << "  [complement flood-fill]\n  H3                        : " << br.H3
      << "  [assumed zero for fallback]\n  Euler chi                 : "
      << br.result_chi << "  [cubical cell count]\n";
  else
    o << "Topology from GUDHI at filtration threshold 0\n  H0 components       "
         ": "
      << br.H0 << "\n  H1 tunnels          : " << br.H1
      << "\n  H2 cavities         : " << br.H2
      << "\n  H3                  : " << br.H3
      << "\n  Euler chi           : " << br.result_chi << "\n";
  o << "\nSurface diagnostics\n  Surface voxels      : " << br.surface_voxels
    << "\n  Surface faces       : " << br.surface_faces
    << "\n\nContact audit\n  H0 face-adjacency / 6-neighbor       : "
    << br.H0_face
    << "\n  H0 closed-cube / 26-neighbor         : " << br.H0_closed
    << "\n  H0 GUDHI                             : ";
  if (br.gudhi_used)
    o << br.H0_gudhi;
  else
    o << "not available; GUDHI disabled";
  o << "\n  Component reduction by contacts      : "
    << br.contact_audit.component_reduction
    << "\n  Inter-component edge contacts        : "
    << br.contact_audit.inter_face_component_edge_contact_count
    << "\n  Inter-component vertex contacts      : "
    << br.contact_audit.inter_face_component_vertex_contact_count
    << "\n  Merger records written               : "
    << br.contact_audit.face_component_merger_count
    << "\n  Contact merger edge records          : "
    << br.contact_audit.merger_edge_record_count
    << "\n  Contact merger vertex records        : "
    << br.contact_audit.merger_vertex_record_count
    << "\n  Contact merger file                  : "
    << contact_audit_file_label(write_contact_files,
                                br.contact_audit.contact_merger_file)
    << "\n  Contact records file                 : "
    << contact_audit_file_label(write_contact_files,
                                br.contact_audit.contact_records_file)
    << "\n\n";
  if (fb)
    o << "Fallback diagnostic sources\n  H0 closed-cube / 26-neighbor         "
         ": "
      << br.H0_closed
      << "\n  H0 face-adjacency / 6-neighbor       : " << br.H0_face
      << "\n  H2 complement                        : " << br.H2_complement
      << "\n  Cell counts C0-C3                    : " << br.C0 << " " << br.C1
      << " " << br.C2 << " " << br.C3
      << "\n  Euler from cells                     : " << br.chi_cells
      << "\n  H1 inference                         : H1 = H0_closed + H2 - chi "
         "= "
      << br.H1
      << "\n\nInternal consistency checks\n  Euler consistency                 "
         "   : "
      << (br.result_chi == br.chi_cells ? "PASS" : "FAIL")
      << "\n  Contact audit consistency            : "
      << (br.errors.empty() ? "PASS" : "FAIL") << "\n";
  else
    o << "Independent cross-checks\n  H0 closed-cube check                 : "
      << (br.H0_closed_check ? "PASS" : "FAIL")
      << "\n  H0 face-adjacency note               : " << h0_face_note(br)
      << "\n  Euler from cells                     : " << br.chi_cells
      << "\n  Euler check                          : "
      << (br.result_chi == br.chi_cells ? "PASS" : "FAIL")
      << "\n  H2 complement                        : " << br.H2_complement
      << "\n  H2 complement check                  : "
      << (br.H2_complement == br.H2 ? "PASS" : "WARN") << "\n";
  for (const auto &e : br.errors)
    o << "Error                 : " << e << "\n";
  o << "\nReady for Stage 6 persistent homology: " << br.readiness_stage6
    << "\n============================================================\n";
  return o.str();
}
static void write_betti_outputs(const std::string &input,
                                const ValidateOptions &vo,
                                const BettiOptions &bo,
                                const ValidationResult &vr,
                                const BettiResult &br, const std::string &rep) {
  std::filesystem::path od(bo.out_dir);
  std::filesystem::create_directories(od);
  bool write_contact_files = should_write_contact_audit_files(bo, br);
  std::vector<std::string> ns = {
      "captop_betti_summary.json", "captop_betti_summary.csv",
      "captop_betti_report.txt", "captop_contact_audit_summary.json"};
  if (write_contact_files) {
    ns.push_back(br.contact_audit.contact_merger_file);
    ns.push_back(br.contact_audit.contact_records_file);
  }
  for (auto &n : ns)
    if (!bo.overwrite && std::filesystem::exists(od / n))
      throw std::runtime_error("output file already exists: " +
                               (od / n).string());
  {
    std::ofstream f(od / "captop_betti_report.txt");
    f << rep;
  }
  if (write_contact_files) {
    std::ofstream m(od / br.contact_audit.contact_merger_file);
    m << "merger_id,face_component_a,face_component_b,closed_component_before_"
         "a,closed_component_before_b,closed_component_after,contact_type,i1,"
         "j1,k1,i2,j2,k2,element_id_1,element_id_2,linear_index_1,linear_index_"
         "2\n";
    for (auto &r : br.contact_audit.merger_records)
      m << r.merger_id << "," << r.face_component_a << "," << r.face_component_b
        << "," << r.closed_component_before_a << ","
        << r.closed_component_before_b << "," << r.closed_component_after << ","
        << contact_type_name(r.contact_type) << "," << r.i1 << "," << r.j1
        << "," << r.k1 << "," << r.i2 << "," << r.j2 << "," << r.k2 << ","
        << r.element_id_1 << "," << r.element_id_2 << "," << r.linear_index_1
        << "," << r.linear_index_2 << "\n";
  }
  if (write_contact_files) {
    std::ofstream cf(od / br.contact_audit.contact_records_file);
    cf << "contact_id,face_component_a,face_component_b,contact_type,i1,j1,k1,"
          "i2,j2,k2,element_id_1,element_id_2,linear_index_1,linear_index_2\n";
    long long id = 1;
    for (auto &r : br.contact_audit.inter_component_contacts)
      cf << id++ << "," << r.component_face_a << "," << r.component_face_b
         << "," << contact_type_name(r.contact_type) << "," << r.i1 << ","
         << r.j1 << "," << r.k1 << "," << r.i2 << "," << r.j2 << "," << r.k2
         << "," << r.element_id_1 << "," << r.element_id_2 << ","
         << r.linear_index_1 << "," << r.linear_index_2 << "\n";
  }
  auto gj = [&](std::ostream &j) {
    j << std::setprecision(17)
      << "{\n  \"computed\": true,\n  \"primary_connectivity\": "
         "\"closed_cubical_complex\",\n  \"face_adjacency_connectivity\": "
         "\"6-neighbor\",\n  \"closed_cube_connectivity\": \"26-neighbor\",\n  "
         "\"H0_face\": "
      << br.H0_face << ",\n  \"H0_closed\": " << br.H0_closed
      << ",\n  \"H0_gudhi\": ";
    if (br.gudhi_used)
      j << br.H0_gudhi;
    else
      j << "null";
    j << ",\n  \"component_reduction\": "
      << br.contact_audit.component_reduction
      << ",\n  \"face_contact_count\": " << br.contact_audit.face_contact_count
      << ",\n  \"edge_contact_count\": " << br.contact_audit.edge_contact_count
      << ",\n  \"vertex_contact_count\": "
      << br.contact_audit.vertex_contact_count
      << ",\n  \"inter_face_component_contact_count\": "
      << br.contact_audit.inter_face_component_contact_count
      << ",\n  \"inter_face_component_edge_contact_count\": "
      << br.contact_audit.inter_face_component_edge_contact_count
      << ",\n  \"inter_face_component_vertex_contact_count\": "
      << br.contact_audit.inter_face_component_vertex_contact_count
      << ",\n  \"inter_face_component_face_contact_count\": "
      << br.contact_audit.inter_face_component_face_contact_count
      << ",\n  \"face_component_merger_count\": "
      << br.contact_audit.face_component_merger_count
      << ",\n  \"merger_edge_record_count\": "
      << br.contact_audit.merger_edge_record_count
      << ",\n  \"merger_vertex_record_count\": "
      << br.contact_audit.merger_vertex_record_count
      << ",\n  \"merger_face_record_count\": "
      << br.contact_audit.merger_face_record_count
      << ",\n  \"gudhi_H0_matches_closed\": ";
    if (br.gudhi_used)
      j << (br.H0_closed_check ? "true" : "false");
    else
      j << "null";
    j << ",\n  \"gudhi_H0_matches_face\": ";
    if (br.gudhi_used)
      j << (br.H0_face_check ? "true" : "false");
    else
      j << "null";
    j << ",\n  \"contact_merger_file\": \""
      << br.contact_audit.contact_merger_file
      << "\",\n  \"contact_records_file\": \""
      << br.contact_audit.contact_records_file
      << "\",\n  \"contact_merger_file_written\": "
      << (write_contact_files ? "true" : "false")
      << ",\n  \"contact_records_file_written\": "
      << (write_contact_files ? "true" : "false") << "\n}";
  };
  {
    std::ofstream aj(od / "captop_contact_audit_summary.json");
    gj(aj);
  }
  {
    std::ofstream j(od / "captop_betti_summary.json");
    j << std::setprecision(17)
      << "{\n  \"software\": {\"name\": \"captop\", \"version\": \""
      << CAPTOP_VERSION << "\"},\n  \"input\": {\"path\": \""
      << json_escape(input) << "\"},\n  \"validation\": {\"grid_mode\": \""
      << grid_mode_name(vo.grid_mode) << "\", \"tolerance\": " << vo.tol
      << ", \"status\": \"VALID\"},\n  \"topology_convention\": "
         "{\"primary_topology\": \"closed_occupied_cubical_complex\", "
         "\"primary_connectivity\": \"26-neighbor closed-cube connectivity\", "
         "\"diagnostic_connectivity\": \"6-neighbor face adjacency\", "
         "\"gudhi_expected_to_match\": \"H0_closed\"},\n  \"topology_engine\": "
         "{\"gudhi_compiled\": "
      << (br.gudhi_compiled ? "true" : "false")
      << ", \"gudhi_used\": " << (br.gudhi_used ? "true" : "false")
      << ", \"gudhi_success\": " << (br.gudhi_success ? "true" : "false")
      << ", \"fallback_used\": true, \"topology_source\": \""
      << topology_source_value(br.topology_source)
      << "\", \"topology_authoritative\": "
      << (br.topology_authoritative ? "true" : "false")
      << ", \"result_authority\": \"" << json_escape(br.result_authority)
      << "\"},\n  \"grid\": {\"nx\": " << vr.nx << ", \"ny\": " << vr.ny
      << ", \"nz\": " << vr.nz << ", \"total_voxels\": " << br.total_voxels
      << ", \"occupied_voxels\": " << br.occupied_voxels
      << ", \"missing_voxels\": " << br.missing_voxels
      << ", \"occupied_fraction\": " << br.occupied_fraction
      << ", \"origin\": [" << vr.origin.x << ", " << vr.origin.y << ", "
      << vr.origin.z << "], \"spacing\": [" << vr.spacing.x << ", "
      << vr.spacing.y << ", " << vr.spacing.z
      << "], \"index_order\": \"i + nx * (j + ny * k)\"},\n  \"topology\": "
         "{\"H0\": "
      << br.H0 << ", \"H1\": " << br.H1 << ", \"H2\": " << br.H2
      << ", \"H3\": " << br.H3 << ", \"chi\": " << br.result_chi
      << ", \"H0_source\": \"" << br.H0_source << "\", \"H1_source\": \""
      << br.H1_source << "\", \"H2_source\": \"" << br.H2_source
      << "\", \"H3_source\": \"" << br.H3_source << "\", \"chi_source\": \""
      << br.chi_source << "\", \"surface_voxels\": " << br.surface_voxels
      << ", \"surface_faces\": " << br.surface_faces
      << "},\n  \"cell_counts\": {\"C0_vertices\": " << br.C0
      << ", \"C1_edges\": " << br.C1 << ", \"C2_faces\": " << br.C2
      << ", \"C3_voxels\": " << br.C3 << ", \"chi_cells\": " << br.chi_cells
      << "},\n  \"cross_checks\": {\"H0_closed\": " << br.H0_closed
      << ", \"H0_face\": " << br.H0_face << ", \"H0_gudhi\": ";
    if (br.gudhi_used)
      j << br.H0_gudhi;
    else
      j << "null";
    j << ", \"H0_closed_check\": \""
      << (br.gudhi_used ? (br.H0_closed_check ? "PASS" : "FAIL") : "NA")
      << "\", \"H0_face_check\": \""
      << (br.gudhi_used ? (br.H0_face_check ? "PASS" : "DIFFERS_EXPECTED")
                        : "DIAGNOSTIC_ONLY")
      << "\", \"chi_result\": " << br.result_chi
      << ", \"chi_matches_cell_count\": "
      << (br.result_chi == br.chi_cells ? "true" : "false")
      << ", \"H2_complement\": " << br.H2_complement
      << ", \"H2_complement_matches_result\": "
      << (br.H2_complement == br.H2 ? "true" : "false")
      << "},\n  \"contact_audit\": ";
    gj(j);
    j << ",\n  \"readiness\": {\"stage6_persistent_homology\": "
      << (br.gudhi_compiled ? "true" : "false") << ", \"reason\": \""
      << (br.gudhi_compiled
              ? "GUDHI enabled"
              : "GUDHI disabled; rebuild with -DCAPTOP_WITH_GUDHI=ON")
      << "\"},\n  \"warnings\": " << json_string_array(br.warnings)
      << ",\n  \"status\": \""
      << (br.success ? (br.gudhi_used ? "VALID AND GUDHI-ANALYZED"
                                      : "VALID AND DIAGNOSTICALLY ANALYZED")
                     : "FAILED")
      << "\"\n}\n";
  }
  {
    std::ofstream c(od / "captop_betti_summary.csv");
    c << "mesh_name,grid_mode,tolerance,field,primary_topology,primary_"
         "connectivity,diagnostic_connectivity,nx,ny,nz,total_voxels,occupied_"
         "voxels,missing_voxels,occupied_fraction,surface_voxels,surface_faces,"
         "H0,H0_source,H0_closed,H0_face,H0_gudhi,H0_closed_check,H0_face_"
         "check,H1,H1_source,H2,H3,chi,C0,C1,C2,C3,chi_cells,component_"
         "reduction,face_contact_count,edge_contact_count,vertex_contact_count,"
         "inter_face_component_contact_count,inter_face_component_edge_contact_"
         "count,inter_face_component_vertex_contact_count,inter_face_component_"
         "face_contact_count,face_component_merger_count,merger_edge_record_"
         "count,merger_vertex_record_count,merger_face_record_count,contact_"
         "merger_file,contact_records_file,contact_merger_file_written,contact_"
         "records_file_written,status,gudhi_"
         "compiled,gudhi_used,gudhi_success,fallback_used,topology_source,"
         "topology_authoritative,result_authority,warnings\n"
      << csv_escape(std::filesystem::path(input).filename().string()) << ","
      << grid_mode_name(vo.grid_mode) << "," << vo.tol << "," << bo.field
      << ",closed_occupied_cubical_complex,"
      << csv_escape("26-neighbor closed-cube connectivity") << ","
      << csv_escape("6-neighbor face adjacency") << "," << vr.nx << "," << vr.ny
      << "," << vr.nz << "," << br.total_voxels << "," << br.occupied_voxels
      << "," << br.missing_voxels << "," << br.occupied_fraction << ","
      << br.surface_voxels << "," << br.surface_faces << "," << br.H0 << ","
      << br.H0_source << "," << br.H0_closed << "," << br.H0_face << ","
      << (br.gudhi_used ? std::to_string(br.H0_gudhi) : "NA") << ","
      << (br.gudhi_used ? (br.H0_closed_check ? "PASS" : "FAIL") : "NA") << ","
      << (br.gudhi_used ? (br.H0_face_check ? "PASS" : "DIFFERS_EXPECTED")
                        : "DIAGNOSTIC_ONLY")
      << "," << br.H1 << "," << br.H1_source << "," << br.H2 << "," << br.H3
      << "," << br.result_chi << "," << br.C0 << "," << br.C1 << "," << br.C2
      << "," << br.C3 << "," << br.chi_cells << ","
      << br.contact_audit.component_reduction << ","
      << br.contact_audit.face_contact_count << ","
      << br.contact_audit.edge_contact_count << ","
      << br.contact_audit.vertex_contact_count << ","
      << br.contact_audit.inter_face_component_contact_count << ","
      << br.contact_audit.inter_face_component_edge_contact_count << ","
      << br.contact_audit.inter_face_component_vertex_contact_count << ","
      << br.contact_audit.inter_face_component_face_contact_count << ","
      << br.contact_audit.face_component_merger_count << ","
      << br.contact_audit.merger_edge_record_count << ","
      << br.contact_audit.merger_vertex_record_count << ","
      << br.contact_audit.merger_face_record_count << ","
      << csv_escape(contact_audit_file_label(
             write_contact_files, br.contact_audit.contact_merger_file))
      << ","
      << csv_escape(contact_audit_file_label(
             write_contact_files, br.contact_audit.contact_records_file))
      << "," << (write_contact_files ? "true" : "false") << ","
      << (write_contact_files ? "true" : "false") << ","
      << csv_escape(br.success
                        ? (br.gudhi_used ? "VALID AND GUDHI-ANALYZED"
                                         : "VALID AND DIAGNOSTICALLY ANALYZED")
                        : "FAILED")
      << "," << (br.gudhi_compiled ? "true" : "false") << ","
      << (br.gudhi_used ? "true" : "false") << ","
      << (br.gudhi_success ? "true" : "false") << ",true,"
      << topology_source_value(br.topology_source) << ","
      << (br.topology_authoritative ? "true" : "false") << ","
      << csv_escape(br.result_authority) << ","
      << csv_escape(br.warnings.empty() ? "" : br.warnings.front()) << "\n";
  }
}
static void print_betti_usage(std::ostream &os) {
  os << "Usage:\n  captop betti <input-mesh> [options]\n\nComputes plain "
        "topology descriptors for binary occupied domains.\n\nIf CAPTOP is "
        "built with GUDHI, Betti numbers are computed from\nGUDHI persistent "
        "Betti numbers at filtration threshold 0.\n\nIf CAPTOP is built "
        "without GUDHI, CAPTOP can produce diagnostic\nfallback estimates:\n  "
        "H0 from closed-cube 26-neighbor union-find\n  H0_face from 6-neighbor "
        "face-adjacency union-find\n  H2 from complement flood-fill\n  chi "
        "from cubical cell counts\n  H1 inferred from H0_closed + H2 - "
        "chi\n\nUse --require-gudhi to reject fallback-only analysis.\n\n";
  print_validation_options(os);
  os << "\nBetti options:\n  --field <prime>          Coefficient field "
        "[default: 2]\n  --out <dir>              Output directory [default: "
        "captop_betti_out]\n  --overwrite              Replace existing Betti "
        "output files\n  --max-memory-gb <value>  Dense-grid memory limit "
        "[default: 2]\n  --force                  Continue above memory "
        "limit\n  --require-gudhi          Fail unless GUDHI support is "
        "compiled and used\n  --allow-fallback         Allow fallback "
        "diagnostics when GUDHI is unavailable [default]\n  --write-diagram    "
        "      Reserve persistence-diagram output (Stage 5 debug)\n  "
        "--write-contact-audit    Write detailed contact audit CSV files\n  "
        "--write-json             Write JSON summary (default outputs include "
        "JSON)\n  --write-csv              Write CSV summary (default outputs "
        "include CSV)\n  --quiet                  Suppress terminal report\n  "
        "-h, --help               Show this help message\n\nContact audit "
        "terms:\n  raw contact              Edge/vertex contact between "
        "distinct 6-neighbor face components\n  merger record            "
        "Independent raw contact that reduces H0_face toward H0_closed\n  "
        "component reduction      H0_face - H0_closed; should equal merger "
        "record count\n";
}
static int run_betti(int argc, char **argv) {
  if (argc < 3) {
    print_betti_usage(std::cerr);
    return 1;
  }
  if (std::string(argv[2]) == "--help" || std::string(argv[2]) == "-h") {
    print_betti_usage(std::cout);
    return 0;
  }
  std::string input = argv[2];
  ParseOptions po;
  ValidateOptions vo;
  BettiOptions bo;
  for (int i = 3; i < argc; ++i) {
    std::string a = argv[i];
    auto need = [&](const std::string &n) {
      if (i + 1 >= argc)
        throw std::runtime_error(n + " requires a value");
      return std::string(argv[++i]);
    };
    try {
      if (a == "--grid") {
        auto v = to_lower(need(a));
        if (v == "strict" || v == "strict-cube")
          vo.grid_mode = GridMode::StrictCube;
        else if (v == "rectilinear")
          vo.grid_mode = GridMode::Rectilinear;
        else if (v == "approximate-cube")
          vo.grid_mode = GridMode::ApproximateCube;
        else
          throw std::runtime_error("unsupported grid mode");
      } else if (a == "--tol") {
        if (!parse_double(need(a), vo.tol) || vo.tol <= 0)
          throw std::runtime_error("invalid --tol");
      } else if (a == "--ignore-non-hexa")
        po.ignore_non_hexa = true;
      else if (a == "--max-errors") {
        long long n;
        if (!parse_long_long(need(a), n) || n <= 0)
          throw std::runtime_error("invalid --max-errors");
        vo.max_errors = n;
      } else if (a == "--field") {
        long long f;
        if (!parse_long_long(need(a), f) || f > INT32_MAX)
          throw std::runtime_error("invalid --field");
        bo.field = f;
      } else if (a == "--out")
        bo.out_dir = need(a);
      else if (a == "--overwrite")
        bo.overwrite = true;
      else if (a == "--force")
        bo.force = true;
      else if (a == "--quiet")
        bo.quiet = true;
      else if (a == "--require-gudhi")
        bo.require_gudhi = true;
      else if (a == "--allow-fallback")
        bo.allow_fallback = true;
      else if (a == "--write-json")
        bo.write_json = true;
      else if (a == "--write-csv")
        bo.write_csv = true;
      else if (a == "--write-diagram")
        bo.write_diagram = true;
      else if (a == "--write-contact-audit")
        bo.write_contact_audit = true;
      else if (a == "--max-memory-gb") {
        if (!parse_double(need(a), bo.max_memory_gb) || bo.max_memory_gb <= 0)
          throw std::runtime_error("invalid --max-memory-gb");
      } else
        throw std::runtime_error("unknown option '" + a + "'");
    } catch (const std::exception &e) {
      std::cerr << "error: " << e.what() << "\n";
      return 1;
    }
  }
  if (!is_prime_field_betti(bo.field)) {
    std::cerr << "error: --field must be a prime integer\n";
    return 1;
  }
#ifndef CAPTOP_WITH_GUDHI
  if (bo.require_gudhi) {
    std::cerr << "CAPTOP was built without GUDHI support; rebuild with "
                 "-DCAPTOP_WITH_GUDHI=ON or omit --require-gudhi to allow "
                 "fallback diagnostics.\n";
    return 1;
  }
#endif
    try{ Mesh mesh=parse_mesh_file(input,po); ValidationResult vr=validate_mesh(mesh,vo); if(!vr.valid){print_validation_report(input,vo,vr); return 2;} BettiResult br=analyze_betti(vr,bo); std::string rep=betti_report(input,vo,bo,vr,br); write_betti_outputs(input,vo,bo,vr,br,rep); if(!bo.quiet) std::cout<<rep; return br.success?0:4; }catch(const ParseError& e){std::cerr<<"parse error: "<<e.what()<<"\n";return 1;}catch(const std::exception& e){std::cerr<<"betti error: "<<e.what()<<"\n";return 3;}
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

    if (command == "convert") {
        return captop::run_convert(argc, argv);
    }

    if (command == "betti") {
        return captop::run_betti(argc, argv);
    }

    if (command == "persist") {
        return captop::run_persist(argc, argv);
    }

    std::cerr << "error: unknown command '" << command << "'\n";
    captop::print_usage(std::cerr);
    return 1;
}
