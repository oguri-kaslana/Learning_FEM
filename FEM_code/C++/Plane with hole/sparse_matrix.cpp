#include "fem_all.h"

// ============================================================
// 1. CooEntry
// ============================================================

CooEntry::CooEntry()
    : row(0), col(0), value(0.0) {}

CooEntry::CooEntry(int row_in, int col_in, double value_in)
    : row(row_in), col(col_in), value(value_in) {}

// ============================================================
// 2. SparseCOO
// ============================================================

SparseCOO::SparseCOO()
    : n_row(0), n_col(0) {}

SparseCOO::SparseCOO(int n_row_in, int n_col_in)
    : n_row(n_row_in), n_col(n_col_in) {}

void SparseCOO::resize_matrix(int n_row_in, int n_col_in) {
    n_row = n_row_in;
    n_col = n_col_in;
    data.clear();
}

void SparseCOO::add_value(int row, int col, double value) {
    // 有限元装配中可能出现非常小的数值。
    // 这里不主动丢弃，避免影响结果精度；转换 CSR 时会自然合并。
    if (row < 0 || row >= n_row || col < 0 || col >= n_col) {
        std::cerr << "[SparseCOO warning] add_value index out of range: row="
                  << row << ", col=" << col << std::endl;
        return;
    }
    data.push_back(CooEntry(row, col, value));
}

void SparseCOO::clear_data() {
    data.clear();
}

void SparseCOO::reserve_data(std::size_t n_data) {
    data.reserve(n_data);
}

int SparseCOO::get_nnz_raw() const {
    return static_cast<int>(data.size());
}

// ============================================================
// 3. SparseCSR
// ============================================================

SparseCSR::SparseCSR()
    : n_row(0), n_col(0) {}

SparseCSR::SparseCSR(int n_row_in, int n_col_in)
    : n_row(n_row_in), n_col(n_col_in) {
    row_ptr.assign(n_row + 1, 0);
}

void SparseCSR::resize_matrix(int n_row_in, int n_col_in) {
    n_row = n_row_in;
    n_col = n_col_in;
    row_ptr.assign(n_row + 1, 0);
    col_id.clear();
    val.clear();
    diag_pos.assign(n_row, -1);
}

void SparseCSR::clear_data() {
    n_row = 0;
    n_col = 0;
    row_ptr.clear();
    col_id.clear();
    val.clear();
    diag_pos.clear();
}

void SparseCSR::multiply_vector(const std::vector<double>& x,
                                std::vector<double>& y) const {
    if (static_cast<int>(x.size()) != n_col) {
        throw std::runtime_error("SparseCSR::multiply_vector: x size does not match matrix column size.");
    }

    y.assign(n_row, 0.0);

    for (int i = 0; i < n_row; ++i) {
        double sum = 0.0;
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
            sum += val[k] * x[col_id[k]];
        }
        y[i] = sum;
    }
}

void SparseCSR::create_diag_pos() {
    diag_pos.assign(n_row, -1);

    for (int i = 0; i < n_row; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
            if (col_id[k] == i) {
                diag_pos[i] = k;
                break;
            }
        }
    }
}

double SparseCSR::get_diag_value(int row) const {
    if (row < 0 || row >= n_row) {
        return 0.0;
    }

    if (row < static_cast<int>(diag_pos.size()) && diag_pos[row] >= 0) {
        return val[diag_pos[row]];
    }

    // 如果没有提前建立 diag_pos，则退化为扫描本行。
    for (int k = row_ptr[row]; k < row_ptr[row + 1]; ++k) {
        if (col_id[k] == row) {
            return val[k];
        }
    }

    return 0.0;
}

double SparseCSR::get_value(int row, int col) const {
    if (row < 0 || row >= n_row || col < 0 || col >= n_col) {
        return 0.0;
    }

    for (int k = row_ptr[row]; k < row_ptr[row + 1]; ++k) {
        if (col_id[k] == col) {
            return val[k];
        }
        if (col_id[k] > col) {
            break;
        }
    }

    return 0.0;
}

int SparseCSR::get_nnz() const {
    return static_cast<int>(val.size());
}

bool SparseCSR::check_valid(std::string& msg) const {
    if (n_row < 0 || n_col < 0) {
        msg = "CSR matrix size is negative.";
        return false;
    }
    if (static_cast<int>(row_ptr.size()) != n_row + 1) {
        msg = "row_ptr size does not equal n_row + 1.";
        return false;
    }
    if (col_id.size() != val.size()) {
        msg = "col_id size does not equal val size.";
        return false;
    }
    if (!row_ptr.empty() && row_ptr.front() != 0) {
        msg = "row_ptr[0] is not zero.";
        return false;
    }
    if (!row_ptr.empty() && row_ptr.back() != static_cast<int>(val.size())) {
        msg = "row_ptr.back() does not equal nnz.";
        return false;
    }
    for (int i = 0; i < n_row; ++i) {
        if (row_ptr[i] > row_ptr[i + 1]) {
            msg = "row_ptr is not nondecreasing.";
            return false;
        }
    }
    for (std::size_t k = 0; k < col_id.size(); ++k) {
        if (col_id[k] < 0 || col_id[k] >= n_col) {
            msg = "col_id contains index out of range.";
            return false;
        }
    }

    msg = "CSR matrix is valid.";
    return true;
}

void SparseCSR::print_info(const std::string& name) const {
    const double rate = calculate_sparse_rate(*this);
    std::cout << "[" << name << "] n_row = " << n_row
              << ", n_col = " << n_col
              << ", nnz = " << get_nnz()
              << ", sparse_rate = " << rate << std::endl;
}

// ============================================================
// 4. COO 转 CSR
// ============================================================

void convert_coo_csr(const SparseCOO& coo, SparseCSR& csr) {
    csr.resize_matrix(coo.n_row, coo.n_col);

    if (coo.n_row <= 0 || coo.n_col <= 0) {
        return;
    }

    // 复制 COO 数据，用于排序和合并。
    std::vector<CooEntry> entries = coo.data;

    // 移除越界项和数值为 NaN/Inf 的项。
    entries.erase(
        std::remove_if(entries.begin(), entries.end(),
            [&](const CooEntry& e) {
                return e.row < 0 || e.row >= coo.n_row ||
                       e.col < 0 || e.col >= coo.n_col ||
                       !std::isfinite(e.value);
            }),
        entries.end());

    // 按行、列排序，便于合并重复项。
    std::sort(entries.begin(), entries.end(),
        [](const CooEntry& a, const CooEntry& b) {
            if (a.row != b.row) return a.row < b.row;
            return a.col < b.col;
        });

    csr.row_ptr.assign(coo.n_row + 1, 0);
    csr.col_id.clear();
    csr.val.clear();

    csr.col_id.reserve(entries.size());
    csr.val.reserve(entries.size());

    const double tiny = 1.0e-30;
    int current_row = 0;
    std::size_t k = 0;

    while (k < entries.size()) {
        const int row = entries[k].row;
        const int col = entries[k].col;
        double sum = 0.0;

        // 将空行的 row_ptr 填到当前非零元数量。
        while (current_row < row) {
            csr.row_ptr[current_row + 1] = static_cast<int>(csr.val.size());
            ++current_row;
        }

        // 合并相同 row-col 的项。
        while (k < entries.size() && entries[k].row == row && entries[k].col == col) {
            sum += entries[k].value;
            ++k;
        }

        // 极小的合并结果可以丢弃，减少无效非零元。
        if (std::fabs(sum) > tiny) {
            csr.col_id.push_back(col);
            csr.val.push_back(sum);
        }
    }

    // 填充剩余行的 row_ptr。
    while (current_row < coo.n_row) {
        csr.row_ptr[current_row + 1] = static_cast<int>(csr.val.size());
        ++current_row;
    }

    csr.create_diag_pos();
}

// ============================================================
// 5. 提取自由自由度子矩阵
// ============================================================

void extract_free_csr(const SparseCSR& K,
                      const std::vector<int>& free_dof,
                      SparseCSR& Kff) {
    const int nf = static_cast<int>(free_dof.size());
    Kff.resize_matrix(nf, nf);

    if (nf == 0) {
        return;
    }

    // old_to_new[旧自由度编号] = 新自由度编号。
    std::vector<int> old_to_new(K.n_row, -1);
    for (int i = 0; i < nf; ++i) {
        const int old_id = free_dof[i];
        if (old_id >= 0 && old_id < K.n_row) {
            old_to_new[old_id] = i;
        }
    }

    std::vector<int> row_count(nf, 0);

    // 第一遍：统计每一行非零元数量。
    for (int i_new = 0; i_new < nf; ++i_new) {
        const int i_old = free_dof[i_new];
        if (i_old < 0 || i_old >= K.n_row) {
            continue;
        }

        for (int k = K.row_ptr[i_old]; k < K.row_ptr[i_old + 1]; ++k) {
            const int j_old = K.col_id[k];
            const int j_new = (j_old >= 0 && j_old < K.n_row) ? old_to_new[j_old] : -1;
            if (j_new >= 0) {
                row_count[i_new]++;
            }
        }
    }

    Kff.row_ptr.assign(nf + 1, 0);
    for (int i = 0; i < nf; ++i) {
        Kff.row_ptr[i + 1] = Kff.row_ptr[i] + row_count[i];
    }

    const int nnz = Kff.row_ptr[nf];
    Kff.col_id.assign(nnz, 0);
    Kff.val.assign(nnz, 0.0);

    std::vector<int> offset = Kff.row_ptr;

    // 第二遍：写入数据。
    for (int i_new = 0; i_new < nf; ++i_new) {
        const int i_old = free_dof[i_new];
        if (i_old < 0 || i_old >= K.n_row) {
            continue;
        }

        for (int k = K.row_ptr[i_old]; k < K.row_ptr[i_old + 1]; ++k) {
            const int j_old = K.col_id[k];
            const int j_new = (j_old >= 0 && j_old < K.n_row) ? old_to_new[j_old] : -1;
            if (j_new >= 0) {
                const int pos = offset[i_new]++;
                Kff.col_id[pos] = j_new;
                Kff.val[pos] = K.val[k];
            }
        }
    }

    // 原 CSR 每行已按列排序，子矩阵通常仍保持排序。
    // 保险起见，这里对每行做一次局部排序，避免后续 get_value 或迭代法依赖列递增时出错。
    for (int i = 0; i < nf; ++i) {
        const int start = Kff.row_ptr[i];
        const int end = Kff.row_ptr[i + 1];
        std::vector<std::pair<int, double> > row_data;
        row_data.reserve(end - start);

        for (int k = start; k < end; ++k) {
            row_data.push_back(std::make_pair(Kff.col_id[k], Kff.val[k]));
        }

        std::sort(row_data.begin(), row_data.end(),
            [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                return a.first < b.first;
            });

        for (int k = start; k < end; ++k) {
            Kff.col_id[k] = row_data[k - start].first;
            Kff.val[k] = row_data[k - start].second;
        }
    }

    Kff.create_diag_pos();
}

void extract_free_vector(const std::vector<double>& F,
                         const std::vector<int>& free_dof,
                         std::vector<double>& Ff) {
    const int nf = static_cast<int>(free_dof.size());
    Ff.assign(nf, 0.0);

    for (int i = 0; i < nf; ++i) {
        const int old_id = free_dof[i];
        if (old_id >= 0 && old_id < static_cast<int>(F.size())) {
            Ff[i] = F[old_id];
        }
    }
}

void recover_full_vector(const std::vector<double>& Uf,
                         const std::vector<int>& free_dof,
                         int n_dof,
                         std::vector<double>& U) {
    U.assign(n_dof, 0.0);

    const int nf = static_cast<int>(free_dof.size());
    for (int i = 0; i < nf; ++i) {
        const int old_id = free_dof[i];
        if (old_id >= 0 && old_id < n_dof && i < static_cast<int>(Uf.size())) {
            U[old_id] = Uf[i];
        }
    }
}

void convert_csr_dense(const SparseCSR& csr,
                       std::vector<double>& A_dense) {
    A_dense.assign(static_cast<std::size_t>(csr.n_row) * static_cast<std::size_t>(csr.n_col), 0.0);

    for (int i = 0; i < csr.n_row; ++i) {
        for (int k = csr.row_ptr[i]; k < csr.row_ptr[i + 1]; ++k) {
            const int j = csr.col_id[k];
            A_dense[static_cast<std::size_t>(i) * csr.n_col + j] = csr.val[k];
        }
    }
}

double calculate_sparse_rate(const SparseCSR& csr) {
    if (csr.n_row <= 0 || csr.n_col <= 0) {
        return 0.0;
    }

    const double total = static_cast<double>(csr.n_row) * static_cast<double>(csr.n_col);
    return static_cast<double>(csr.get_nnz()) / total;
}
