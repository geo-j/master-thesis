#include <vector>

template<typename type>
void push_back_unique(std::vector<type> &v, type element) {
    auto it = std::find(v.begin(), v.end(), element);

    if (it == v.end())
        v.push_back(element);
}