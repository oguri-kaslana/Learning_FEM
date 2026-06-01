#include "fem_all.h"

using namespace std;

// ============================================================
// 构造函数
// ============================================================
TimeRecord::TimeRecord() {
    time_table.clear();
    start_table.clear();
}

// ============================================================
// 开始计时
// ============================================================
void TimeRecord::start_time(const string& name) {
    start_table[name] = chrono::high_resolution_clock::now();
}

// ============================================================
// 停止计时并累加
// ============================================================
void TimeRecord::stop_time(const string& name) {
    auto it = start_table.find(name);

    if (it == start_table.end()) {
        cerr << "[TimeRecord warning] step \"" << name
             << "\" has not been started, stop_time is ignored." << endl;
        return;
    }

    chrono::high_resolution_clock::time_point end_time =
        chrono::high_resolution_clock::now();

    double seconds = chrono::duration<double>(end_time - it->second).count();

    time_table[name] += seconds;

    // 删除本次开始时间，避免重复 stop 同一步骤产生错误计时
    start_table.erase(it);
}

// ============================================================
// 手动累加时间
// ============================================================
void TimeRecord::add_time(const string& name, double seconds) {
    if (seconds < 0.0) {
        cerr << "[TimeRecord warning] negative time is ignored for step \""
             << name << "\"." << endl;
        return;
    }

    time_table[name] += seconds;
}

// ============================================================
// 获取某一步骤累计耗时
// ============================================================
double TimeRecord::get_time(const string& name) const {
    auto it = time_table.find(name);

    if (it == time_table.end()) {
        return 0.0;
    }

    return it->second;
}

// ============================================================
// 判断某一步骤是否存在计时记录
// ============================================================
bool TimeRecord::has_time(const string& name) const {
    return time_table.find(name) != time_table.end();
}

// ============================================================
// 清空计时器
// ============================================================
void TimeRecord::clear_time() {
    time_table.clear();
    start_table.clear();
}

// ============================================================
// 打印全部计时结果
// ============================================================
void TimeRecord::print_time() const {
    cout << "\n================ Time Summary ================\n";
    cout << left << setw(40) << "step_name"
         << right << setw(18) << "time_seconds" << "\n";
    cout << "----------------------------------------------------------\n";

    for (auto it = time_table.begin(); it != time_table.end(); ++it) {
        cout << left << setw(40) << it->first
             << right << setw(18) << fixed << setprecision(8)
             << it->second << "\n";
    }

    if (!start_table.empty()) {
        cout << "----------------------------------------------------------\n";
        cout << "[TimeRecord warning] Some timers are still running:\n";
        for (auto it = start_table.begin(); it != start_table.end(); ++it) {
            cout << "  - " << it->first << "\n";
        }
    }

    cout << "================================================\n\n";
}

// ============================================================
// 输出计时表到 csv
// ============================================================
void TimeRecord::write_time(const string& filename, const string& mode) const {
    ofstream fout(filename.c_str());

    if (!fout.is_open()) {
        cerr << "[TimeRecord warning] cannot open time file: "
             << filename << endl;
        return;
    }

    fout << "step_name,time_seconds,mode\n";
    fout << fixed << setprecision(12);

    for (auto it = time_table.begin(); it != time_table.end(); ++it) {
        fout << it->first << ","
             << it->second << ","
             << mode << "\n";
    }

    // 如果存在未停止计时器，也写入文件，方便调试。
    // 时间记为 -1，表示该步骤 start 之后没有正常 stop。
    for (auto it = start_table.begin(); it != start_table.end(); ++it) {
        fout << it->first << ","
             << -1.0 << ","
             << mode << "_not_stopped\n";
    }

    fout.close();
}

// ============================================================
// 获取内部计时表
// ============================================================
const map<string, double>& TimeRecord::get_table() const {
    return time_table;
}
