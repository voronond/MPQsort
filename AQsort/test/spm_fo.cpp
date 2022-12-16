#include <algorithm>
#include <iostream>
#include <vector>

#include <aqsort.h>

template <typename T>
struct Comp
{
    Comp(std::vector<T>& rows, std::vector<T>& cols) : rows_(rows), cols_(cols) { }

    inline bool operator()(std::size_t i, std::size_t j) const
    {
        if (rows_[i] < rows_[j])
            return true;
        if ((rows_[i] == rows_[j]) && (cols_[i] < cols_[j]))
            return true;
        return false;
    }

    private: std::vector<T> &rows_, &cols_;
};

template <typename T, typename U>
struct Swap
{
    Swap(std::vector<T>& rows, std::vector<T>& cols, std::vector<U>& vals)
        : rows_(rows), cols_(cols), vals_(vals) { }

    inline void operator()(std::size_t i, std::size_t j)
    {
        std::swap(rows_[i], rows_[j]);
        std::swap(cols_[i], cols_[j]);
        std::swap(vals_[i], vals_[j]);
    }

    private: 
    std::vector<T> &rows_, &cols_;
    std::vector<U> &vals_;
};

int main()
{
    typedef unsigned int uint32_t;

    std::vector<uint32_t> rows, cols;
    std::vector<double>   vals;

    rows.push_back(0); cols.push_back(0); vals.push_back(1.0);
    rows.push_back(3); cols.push_back(0); vals.push_back(5.0);
    rows.push_back(1); cols.push_back(1); vals.push_back(3.0);
    rows.push_back(2); cols.push_back(2); vals.push_back(4.0);
    rows.push_back(0); cols.push_back(3); vals.push_back(2.0);
    rows.push_back(3); cols.push_back(3); vals.push_back(6.0);

    for (int i = 0; i < rows.size(); i++) 
        std::cout << "(" << rows[i] << ", " << cols[i] << ") = " << vals[i] << std::endl;
    std::cout << std::endl;

    Comp<uint32_t>         comp(rows, cols);
    Swap<uint32_t, double> swap(rows, cols, vals);

    aqsort::sort(rows.size(), &comp, &swap);

    for (int i = 0; i < rows.size(); i++) 
        std::cout << "(" << rows[i] << ", " << cols[i] << ") = " << vals[i] << std::endl;
}
