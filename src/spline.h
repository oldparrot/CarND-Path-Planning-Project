#pragma once

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>

namespace {

namespace tk {

// 带状矩阵 solver
class band_matrix {
 public:
    band_matrix() {};

    band_matrix(int dim, int n_u, int n_l);

    ~band_matrix() {};

    void resize(int dim, int n_u, int n_l);

    int dim() const;

    int num_upper() const {
        return m_upper.size() - 1;
    }

    int num_lower() const {
        return m_lower.size() - 1;
    }

    double &operator() (int i, int j);

    double operator() (int i, int j) const;

    double &saved_diag(int i);

    double saved_diag(int i) const;

    void lu_decompose();

    std::vector<double> r_solve(const std::vector<double> &b) const;

    std::vector<double> l_solve(const std::vector<double> &b) const;

    std::vector<double> lu_solve(const std::vector<double> &b,
                                 bool is_lu_decomposed = false);

 private:
    std::vector<std::vector<double>> m_upper;
    std::vector<std::vector<double>> m_lower;
};

class spline {
 public:
    enum bd_type {
        first_deriv = 1,
        second_deriv = 2
    };

 public:
    spline() : m_left(second_deriv), m_right(second_deriv),
               m_left_value(0.0), m_right_value(0.0),
               m_force_linear_extrapolation(false) {}
    
    void set_boundary(bd_type left, double left_value,
                      bd_type right, double right_value,
                      bool force_linear_extrapolation = false);
    
    void set_points(const std::vector<double> &x,
                    const std::vector<double> &y, bool cubic_spline = true);
    
    double operator() (double x) const;

 private:
    std::vector<double> m_x;
    std::vector<double> m_y;
    // f(x) = a * (x - x_i) ^ 3 + b * (x - x_i) ^ 2 + c * (x - x_i) + y_i
    std::vector<double> m_a; // spline coefficients
    std::vector<double> m_b; // spline coefficients
    std::vector<double> m_c; // spline coefficients

    double m_b0;
    double m_c0;

    bd_type m_left;
    bd_type m_right;

    double m_left_value;
    double m_right_value;

    bool m_force_linear_extrapolation;
};


band_matrix::band_matrix(int dim, int n_u, int n_l) {
    resize(dim, n_u, n_l);
}

void band_matrix::resize(int dim, int n_u, int n_l) {
    assert(dim > 0);
    assert(n_u >= 0);
    assert(n_l >= 0);
    m_upper.resize(n_u + 1);
    m_lower.resize(n_l + 1);

    for(size_t i = 0; i < m_upper.size(); i++) {
        m_upper[i].resize(dim);
    }

    for(size_t i = 0; i < m_lower.size(); i++) {
        m_lower[i].resize(dim);
    }
}

int band_matrix::dim() const {
    if(m_upper.size() > 0) {
        return m_upper[0].size();
    } else {
        return 0;
    }
}

double &band_matrix::operator() (int i, int j) {
    int k = j - i;
    assert((i >= 0) && (i < dim()) && (j >= 0) && (j < dim()));
    assert((-num_lower() <= k) && (k <= num_upper()));
    if(k >= 0) {
        return m_upper[k][i];
    } else {
        return m_lower[-k][i];
    }
}

double band_matrix::operator() (int i, int j) const {
    int k = j - i;
    assert((i >= 0) && (i < dim()) && (j >= 0) && (j < dim()));
    assert((-num_lower() <= k) && (k <= num_upper()));
    if(k >= 0) {
        return m_upper[k][i];
    } else {
        return m_lower[-k][i];
    }
}

double &band_matrix::saved_diag(int i) {
    assert((i >= 0) && (i < dim()));
    return m_lower[0][i];
}

double band_matrix::saved_diag(int i) const {
    assert((i >= 0) && (i < dim()));
    return m_lower[0][i];
}

// 带状矩阵的LU分解(可尝试改为LUP分解), 求解后分别存入m_lower, m_upper
void band_matrix::lu_decompose() {
    int i_max, j_max;
    int j_min;
    double x;
    // 将矩阵的对角元素更新为1, m_lower[0]存计算元素， m_upper[0]存矩阵元素
    for(int i = 0; i < this->dim(); i++) {
        assert(this->operator()(i, i) != 0.0);
        this->saved_diag(i) = 1.0 / this->operator()(i,i);
        j_min = std::max(0, i - this->num_lower());
        j_max = std::min(this->dim() - 1, i + this->num_upper());
        for(int j = j_min; j <= j_max; j++) {
            this->operator()(i, j) *= this->saved_diag(i);
        }
        this->operator()(i, i) = 1.0;
    }

    for(int k = 0; k < this->dim(); k++) {
        i_max = std::min(this->dim() - 1, k + this->num_lower());
        for(int i = k + 1; i <= i_max; i++) {
            assert(this->operator()(k, k) != 0.0);
            x = -this->operator()(i, k) / this->operator()(k, k);
            this->operator()(i, k) = -x;
            j_max = std::min(this->dim() - 1, k + this->num_upper());
            for(int j = k + 1; j <= j_max; j++) {
                this->operator()(i, j) = this->operator()(i, j) + x * this->operator()(k, j);
            }
        }
    }
}

// 求解Ly = b
std::vector<double> band_matrix::l_solve(const std::vector<double> &b) const {
    assert(this->dim() == (int)b.size());
    std::vector<double> x(this->dim());
    int j_start;
    double sum;
    for(int i = 0; i < this->dim(); i++) {
        sum = 0;
        j_start = std::max(0, i - this->num_lower());
        for(int j = j_start; j < i; j++) {
            sum += this->operator()(i, j) * x[j];
        }
        x[i] = (b[i] * this->saved_diag(i)) - sum;
    }
    return x;
}

// 求解Ux = y
std::vector<double> band_matrix::r_solve(const std::vector<double> &b) const {
    assert(this->dim() == (int)b.size());
    std::vector<double> x(this->dim());
    int j_stop;
    double sum;
    for(int i = this->dim() - 1; i >= 0; i--) {
        sum = 0;
        j_stop = std::min(this->dim() - 1, i + this->num_upper());
        for(int j = i + 1; j <= j_stop; j++) {
            sum += this->operator()(i, j) * x[j];
        }
        x[i] = (b[i] - sum) / this->operator()(i, i);
    }
    return x;
}

std::vector<double> band_matrix::lu_solve(const std::vector<double> &b,
                                          bool is_lu_decomposed) {
    assert(this->dim() == (int)b.size());
    std::vector<double> x;
    std::vector<double> y;
    if(!is_lu_decomposed) {
        this->lu_decompose();
    }
    y = this->l_solve(b);
    x = this->r_solve(y);
    return x;
}

// spline implementation
void spline::set_boundary(spline::bd_type left, double left_value,
                          spline::bd_type right, double right_value,
                          bool force_linear_extrapolation) {
    m_left = left;
    m_right = right;
    m_left_value = left_value;
    m_right_value = right_value;
    m_force_linear_extrapolation = force_linear_extrapolation;
}

/* 
https://zhuanlan.zhihu.com/p/62860859 
不同于知乎的地方在于这里m_i = c_i， 因此系数为1/3, 2/3, 1/3
*/
void spline::set_points(const std::vector<double> &x,
                        const std::vector<double> &y, bool cubic_spline) {
    assert(x.size() == y.size());
    assert(x.size() > 2);
    m_x = x;
    m_y = y;
    int n = x.size();
    for(int i = 0; i < n - 1; i++) {
        assert(m_x[i] < m_x[i + 1]);
    }

    if(cubic_spline) {
        band_matrix A(n, 1, 1);
        std::vector<double> rhs(n);
        for(int i = 1; i < n - 1; i++) {
            A(i, i - 1) = 1.0 / 3.0 * (x[i] - x[i-1]);
            A(i, i) = 2.0 / 3.0 *(x[i + 1] - x[i - 1]);
            A(i, i + 1) = 1.0 / 3.0 *(x[i + 1] - x[i]);
            rhs[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) 
                   - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        }

        if(m_left == spline::second_deriv) {
            A(0, 0) = 2.0;
            A(0, 1) = 0.0;
            rhs[0] = m_left_value;
        } else if(m_left == spline::first_deriv) {
            A(0, 0) = 2.0 * (x[1] - x[0]);
            A(0, 1) = 1.0 * (x[1] - x[0]);
            rhs[0] = 3.0 * ((y[1] - y[0]) / (x[1] - x[0]) - m_left_value);
        } else {
            assert(false);
        }

        if(m_right == spline::second_deriv) {
            // 2*b[n-1] = f''
            A(n-1,n-1)=2.0;
            A(n-1,n-2)=0.0;
            rhs[n-1]=m_right_value;
        } else if(m_right == spline::first_deriv) {
            // c[n-1] = f', needs to be re-expressed in terms of b:
            // (b[n-2]+2b[n-1])(x[n-1]-x[n-2])
            // = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
            A(n-1,n-1)=2.0*(x[n-1]-x[n-2]);
            A(n-1,n-2)=1.0*(x[n-1]-x[n-2]);
            rhs[n-1]=3.0*(m_right_value-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
        } else {
            assert(false);
        }

        // solve the equation system to obtain the parameters b[]
        m_b=A.lu_solve(rhs);

        // calculate parameters a[] and c[] based on b[]
        m_a.resize(n);
        m_c.resize(n);
        for(int i=0; i<n-1; i++) {
            m_a[i]=1.0/3.0*(m_b[i+1]-m_b[i])/(x[i+1]-x[i]);
            m_c[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
                   - 1.0/3.0*(2.0*m_b[i]+m_b[i+1])*(x[i+1]-x[i]);
        }
    } else { // linear interpolation
        m_a.resize(n);
        m_b.resize(n);
        m_c.resize(n);
        for(int i=0; i<n-1; i++) {
            m_a[i]=0.0;
            m_b[i]=0.0;
            m_c[i]=(m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i]);
        }
    }

    // for left extrapolation coefficients
    m_b0 = (m_force_linear_extrapolation==false) ? m_b[0] : 0.0;
    m_c0 = m_c[0];

    // for the right extrapolation coefficients
    // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
    double h=x[n-1]-x[n-2];
    // m_b[n-1] is determined by the boundary condition
    m_a[n-1]=0.0;
    m_c[n-1]=3.0*m_a[n-2]*h*h+2.0*m_b[n-2]*h+m_c[n-2];   // = f'_{n-2}(x_{n-1})
    if(m_force_linear_extrapolation==true)
        m_b[n-1]=0.0;
}

double spline::operator() (double x) const {
    size_t n = m_x.size();
    std::vector<double>::const_iterator it;
    it = std::lower_bound(m_x.begin(), m_x.end(), x);
    int idx = std::max(int(it - m_x.begin()) - 1, 0);

    double h = x - m_x[idx];
    double interpol;

    if(x < m_x[0]) {
        // extrapolation to the left
        interpol = (m_b0 * h + m_c0) * h + m_y[0];
    } else if(x > m_x[n-1]) {
        // extrapolation to the right
        interpol=(m_b[n-1] * h + m_c[n-1]) * h + m_y[n-1];
    } else {
        // interpolation
        interpol=((m_a[idx] * h + m_b[idx]) * h + m_c[idx]) * h + m_y[idx];
    }
    return interpol;
}

} // namespace tk
} // namespace