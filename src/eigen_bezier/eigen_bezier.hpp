#include <fstream>
#include <vector>
#include <type_traits>
#include <Eigen/Eigen>

namespace eigen_bezier
{
template <typename Integer>
std::enable_if_t<std::is_integral_v<Integer>, Integer> Factorial(Integer n)
{
    if(n == 0)
    {
        return 1;
    }
    else
    {
        Integer result = 1;
        for(Integer i = 2; i <= n; i++)
        {
            result *= i;
        }
        return result;
    }
}

template <typename Scalar>
class BernsteinPolynomial
{
public:
    BernsteinPolynomial(int degree);
    int Coeff(int i) const;
    Scalar Term(int i, Scalar t) const;
    Scalar Eval(Scalar t) const;
    int Degree() const;
protected:
    int BinomialCoeff(int n, int v) const;
    int _degree = 0;
    std::vector<int> _coeffs;
};

template <typename Scalar>
BernsteinPolynomial<Scalar>::BernsteinPolynomial(int degree)
    :_degree(degree)
{
    for(int i = 0; i <= _degree; i++)
    {
        _coeffs.push_back(BinomialCoeff(_degree, i));
    }
}

template <typename Scalar>
int BernsteinPolynomial<Scalar>::Coeff(int i) const
{
    if(i < 0 || i >= _coeffs.size())
    {
        return 0;
    }
    return _coeffs[i];
}

template <typename Scalar>
Scalar BernsteinPolynomial<Scalar>::Term(int i, Scalar t) const
{
    Scalar tn = 1.0;
    Scalar tn_1 = 1.0;
    for(int j = 0; j < i; j++)
    {
        tn *= t;
        
    }
    for(int j = 0; j < (_degree - i); j++)
    {
        tn_1 *= (1.0 - t);
    }
    
    return  tn * tn_1 * _coeffs[i];
}

template <typename Scalar>
Scalar BernsteinPolynomial<Scalar>::Eval(Scalar t) const
{
    t = std::min(std::max(t, 0.0), 1.0);
    Scalar result = 0.0;
    for(int i = 0; i <= _degree; i++)
    {
        result += Term(i, t);
    }
    return result;
}

template <typename Scalar>
int BernsteinPolynomial<Scalar>::Degree() const
{
    return _degree;
}

template <typename Scalar>
int BernsteinPolynomial<Scalar>::BinomialCoeff(int n, int v) const
{
    return Factorial(n) / (Factorial(v) * (Factorial(n - v)));
}


template <typename Scalar>
class Curve
{
public:
    Curve( const std::vector<Eigen::Vector3<Scalar>>& points);
    Eigen::Vector3<Scalar> Eval(float t) const;
    int Degree() const;
    std::vector<Eigen::Vector3<Scalar>> Sample(uint32_t nb_sample) const;
    bool WriteOBJ(const std::string& path) const;
protected:
    std::vector<Eigen::Vector3<Scalar>> _points;
    BernsteinPolynomial<Scalar> _polynomial;
};

template <typename Scalar>
Curve<Scalar>::Curve( const std::vector<Eigen::Vector3<Scalar>>& points )
    :_points(points), _polynomial(points.size() - 1)
{
}

template <typename Scalar>
Eigen::Vector3<Scalar> Curve<Scalar>::Eval(float t) const
{
    Eigen::Vector3<Scalar> Result = Eigen::Vector3<Scalar>::Zero();
    for(int i = 0; i < _points.size(); i++)
    {
        Result += _polynomial.Term(i, t) * _points[i];
    }
    return Result;
}

template <typename Scalar>
int Curve<Scalar>::Degree() const
{
    return _points.size() - 1;
}

template <typename Scalar>
std::vector<Eigen::Vector3<Scalar>> Curve<Scalar>::Sample(uint32_t nb_sample) const
{
    if(nb_sample <= 1)
    {
        return {};
    }
    float step = 1.0 / nb_sample;
    float t = 0.0;
    std::vector<Eigen::Vector3<Scalar>> samples;
    for(uint32_t i = 0; i < nb_sample; i++)
    {
        samples.push_back(Eval(t));
        t += step;
    }
    return samples;
}

template <typename Scalar>
bool Curve<Scalar>::WriteOBJ( const std::string& path ) const
{
    std::vector<Eigen::Vector3<Scalar>> samples = Sample(1000);
    if(samples.empty())
    {
        return false;
    }
    std::ofstream ofs(path);
    if(ofs.fail())
    {
        return false;
    }
    for(size_t i = 0; i < samples.size() - 1; i++)
    {
        ofs << "v " << samples[i][0] << ' ' << samples[i][1] << ' ' << samples[i][2] << "\n";
        ofs << "v " << samples[i + 1][0] << ' ' << samples[i + 1][1] << ' ' << samples[i + 1][2] << "\n";
        ofs << "l " << i * 2 + 1 << ' ' << i * 2 + 2 << "\n";
    }
    ofs.close();
    return ofs.good();
}

}
