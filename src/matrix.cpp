#include <iostream>
#include <string>
#include <random>
#include <boost/numeric/odeint.hpp>

#include "matrix_type.hpp"

const size_t N = 8;
const size_t K = 3;

typedef matrix<double, N> matrix_type;
typedef foo<double, N> f_type;
typedef boost::array<matrix_type, K> container_type;

struct bfss_coor
{
  bfss_coor(void) {}

  void operator()(const container_type &p, container_type &dqdt) const
  {
    for (size_t i = 0; i < K; ++i)
    {
      dqdt[i] = p[i];
    }
  }
};

struct bfss_momentum
{

  const f_type &m_F;

  bfss_momentum(const f_type &F) : m_F(F) {}

  void operator()(const container_type &q, container_type &dpdt) const
  {
    for (size_t k = 0; k < K; ++k)
    {
      dpdt[k] = 0.0;
      for (size_t i = 0; i < K; ++i)
      {
        dpdt[k] += m_F.comm(m_F.comm(q[k], q[i]), q[i]);
      }
    }
  }
};

double radius(const container_type &q)
{
  double sum = 0.0;
  for (size_t i = 0; i < K; ++i)
    sum += 0.5 * norm(q[i]);
  return sum;
}

double energy(const container_type &q, const container_type &p, const f_type &F)
{
  double sum = 0.0;
  for (size_t i = 0; i < K; ++i)
    sum += 0.5 * norm(p[i]);

  for (size_t i = 0; i < K; ++i)
    for (size_t j = i + 1; j < K; ++j)
      sum += 0.5 * norm(F.comm(q[i], q[j]));

  return sum;
}

struct streaming_observer_elements
{
  std::ostream &m_out;

  streaming_observer_elements(std::ostream &out) : m_out(out) {}

  template <class State>
  void operator()(const State &s, double t) const
  {
    container_type q = s.first;
    m_out << t;
    for (size_t i = 0; i < K; ++i)
      m_out << "    " << q[i];
    m_out << "\n";
  }
};

struct streaming_observer_elements_with_momenta
{
  std::ostream &m_out;

  streaming_observer_elements_with_momenta(std::ostream &out) : m_out(out) {}

  template <class State>
  void operator()(const State &s, double t) const
  {
    container_type q = s.first;
    container_type p = s.second;
    m_out << t;
    for (size_t i = 0; i < K; ++i)
      m_out << "    " << q[i];
    m_out << "    ";
    for (size_t i = 0; i < K; ++i)
      m_out << "    " << p[i];
    m_out << "\n";
  }
};

struct streaming_observer_energy
{
  std::ostream &m_out;
  const f_type &m_F;

  streaming_observer_energy(std::ostream &out, const f_type &F) : m_out(out), m_F(F) {}

  template <class State>
  void operator()(const State &s, double t) const
  {
    container_type q = s.first;
    container_type p = s.second;
    m_out << t << " " << energy(q, p, m_F);
    m_out << "\n";
  }
};

struct streaming_observer_radius
{
  std::ostream &m_out;

  streaming_observer_radius(std::ostream &out) : m_out(out) {}

  template <class State>
  void operator()(const State &s, double t) const
  {
    container_type q = s.first;
    m_out << t << " " << radius(q);
    m_out << "\n";
  }
};


boost::array<double, N*N-1> randomState()
{
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<> dist(-1.0, 1.0);

  boost::array<double, N*N-1> state;

  for (size_t i = 0; i < state.size(); ++i)
  {
    state[i] = dist(gen);
  }

  return state;
}


int main(int argc, char **argv)
{
  using namespace std;
  using namespace boost::numeric::odeint;

  const string f_file = "/home/rishi/Projects/bfss-cpp/F/N" + to_string(N) + ".dat";
  const f_type masterF(f_file);

  container_type q = {{
    matrix_type(randomState()),
    matrix_type(randomState()),
    matrix_type(randomState())
  }};

  container_type p = {{
      matrix_type(0.0),
      matrix_type(0.0),
      matrix_type(0.0)
  }};

  typedef symplectic_rkn_sb3a_mclachlan<container_type> stepper_type;

  const double start_time = 0.0;
  const double end_time = 1000.0;
  const double dt = 0.1;

  integrate_const(
      stepper_type(),
      make_pair(bfss_coor(), bfss_momentum(masterF)),
      make_pair(boost::ref(q), boost::ref(p)),
      start_time,
      end_time,
      dt,
      streaming_observer_elements_with_momenta(cout));

  return 0;
}
