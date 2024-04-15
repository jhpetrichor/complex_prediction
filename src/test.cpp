#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

// 自定义比较函数
bool comparePairs(const pair<int, int>& pair1, const pair<int, int>& pair2) {
    return pair1.second < pair2.second;
}

int main() {
    // 假设有一个 vector<pair<Edge*, int>> 对象
    vector<pair<int, int>> edges;

    // 添加一些元素用于示例
    edges.push_back(make_pair(0, 5));
    edges.push_back(make_pair(2, 2));
    edges.push_back(make_pair(10, 8));
    edges.push_back(make_pair(1, 3));

    // 使用自定义的比较函数对 vector 进行排序
    sort(edges.begin(), edges.end(), comparePairs);

    // 输出排序后的结果
    for (const auto& pair : edges) {
        cout << pair.second << " ";
    }
    cout << endl;

    return 0;
}