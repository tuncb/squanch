#pragma once
#include <set>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <quadrature/Gaussian.h>
#include <squanch/Nurbs.h>
#include <squanch/IgaAlgorithms.h>
#include <squanch/PointInversion.h>

#include <mkl.h>

namespace squanch
{
  namespace detail {

		template <typename T> std::vector<T> findSmallestKnotIntervals(const std::vector<T>& knots1, const std::vector<T>& knots2)
		{
			std::set<T> delta_set;
			delta_set.insert(knots1.begin(), knots1.end());
			delta_set.insert(knots2.begin(), knots2.end());
			
			return std::vector<T>(delta_set.begin(), delta_set.end());
		}

    template <template<class, unsigned int> class NurbsType, typename T>
    std::vector<unsigned int> findGlobalPositions(const NurbsType<T, 1>& nurbs, const typename NurbsType<T, 1>::SpanType& span)
    {
      std::vector<unsigned int> positions;
      auto p = nurbs.curves()[0].p();

      for (unsigned int i = 0; i <= p; ++i) {
        positions.push_back(span - p + i);
      }
      return positions;
    }

    template <template<class, unsigned int> class NurbsType, typename T>
    std::vector<unsigned int> findGlobalPositions(const NurbsType<T, 2>& nurbs, const typename NurbsType<T, 2>::SpanType& span)
    {
      std::vector<unsigned int> positions;
      auto p1 = nurbs.curves()[0].p();
      auto p2 = nurbs.curves()[1].p();

      auto&& shape = nurbs.cpi().shape();


      for (unsigned int i = 0; i <= p2; ++i) {
        auto numI = span[1] - p2 + i;
        for (unsigned int j = 0; j <= p1; ++j) {
          positions.push_back(numI * shape[0] + span[0] - p1 + j);
        }
      }
      return positions;
    }

    template <template<class, unsigned int> class NurbsType, typename T>
    std::vector<unsigned int> findGlobalPositions(const NurbsType<T, 3>& nurbs, const typename NurbsType<T, 3>::SpanType& span)
    {
      std::vector<unsigned int> positions;
      auto p1 = nurbs.curves()[0].p();
      auto p2 = nurbs.curves()[1].p();
      auto p3 = nurbs.curves()[2].p();

      auto&& shape = nurbs.cpi().shape();

      for (unsigned int k = 0; k <= p3; ++k) {
        auto numK = span[2] - p3 + k;
        for (unsigned int i = 0; i <= p2; ++i) {
          auto numI = span[1] - p2 + i;
          for (unsigned int j = 0; j <= p1; ++j) {
            positions.push_back(numK * (shape[0] * shape[1]) + numI * shape[0] + span[0] - p1 + j);
          }
        }
      }
      return positions;
    }

    template <typename T>
    void assemble_mortar_matrix(const std::vector<unsigned int>& positionsH, const std::vector<unsigned int>& positionsG, const Eigen::Matrix<T, -1, -1>& local_m, Eigen::Matrix<T, -1, -1>& global_m)
    {
      for (auto i = 0u; i < positionsG.size(); ++i) {
        auto x = positionsG[i];
        for (auto j = 0u; j < positionsH.size(); ++j) {
          auto y = positionsH[j];
          global_m(x, y) += local_m(i, j);
        }
      }
    }

    struct IdentityPointReverser
    {
      template <template<class, unsigned int> class NurbsType, typename T, unsigned int NDIM>
      void operator()(const squanch::Model<T>&, const NurbsType<T, NDIM>&, const NurbsType<T, NDIM>&,
        typename const NurbsType<T, NDIM>::SpanType& spanH, typename const NurbsType<T, NDIM>::ParametricPointType& uH,
        typename NurbsType<T, NDIM>::SpanType& spanG, typename NurbsType<T, NDIM>::ParametricPointType& uG)
      {
        spanG = spanH;
        uG = uH;
      }
    };


		inline void cout_pp(double h) { std::cout << h << "\n"; }
		inline void cout_pp(std::array<double, 2> h) { std::cout << h[0] << "," << h[1] << "\n"; }



    struct NewtonRaphsonPointReverser
    {
      template <template<class, unsigned int> class NurbsType, typename T, unsigned int DimH, unsigned int DimG>
      void operator()(const squanch::Model<T>& model, const NurbsType<T, DimH>& host, const NurbsType<T, DimG>& guest,
        typename const NurbsType<T, DimH>::SpanType& spanH, typename const NurbsType<T, DimH>::ParametricPointType& uH,
        typename NurbsType<T, DimG>::SpanType& spanG, typename NurbsType<T, DimG>::ParametricPointType& uG)
      {
        // find guest point
        squanch::PointContainerTypes<T>::point_type P;
        Eigen::Matrix<T, -1, 1> weights(num_nodes_per_span(host.curves()));
        Eigen::Matrix<T, -1, 1> r(num_nodes_per_span(host.curves()));
        gather_weights(model, host, spanH, weights);
        compute_b_basis(host, spanH, uH, r);
        P = compute_nurbs_point(model, host, uH, spanH, r);
        // reverse host point
        uG = inverse_point_nr(model, guest, uH, P);
				spanG = find_span(guest, uG);
      }
    };


    template <template<class, unsigned int> class NurbsType, typename T, unsigned int DimG, typename PointReverserType>
    void mortar_integrator(const squanch::Model<T>& model, PointReverserType& pointReverser, const ozp::quadrature::Quadrature& quad, const NurbsType<T, 1>& host, const NurbsType<T, DimG>& guest, 
			Eigen::Matrix<T, -1, -1>& mortar)
    {
      const T small_number = 0.0000001;

			mortar.setZero();

      auto num_local_nodes_host = num_nodes_per_span(host.curves());
      auto num_local_nodes_guest = num_nodes_per_span(guest.curves());

      Eigen::Matrix<T, -1, -1> mortar_l(num_local_nodes_guest, num_local_nodes_host);
      Eigen::Matrix<T, -1, 1> host_weights(num_local_nodes_host);
      Eigen::Matrix<T, -1, 1> guest_weights(num_local_nodes_guest);
      Eigen::Matrix<T, -1, 1> host_r(num_local_nodes_host);
      Eigen::Matrix<T, -1, 1> guest_r(num_local_nodes_guest);

      NurbsType<T, DimG>::ParametricPointType uG;
      NurbsType<T, DimG>::SpanType spanG;

			auto knotIntervals = detail::findSmallestKnotIntervals(host.curves()[0].knots(), guest.curves()[0].knots());

			for (auto i = 1u; i < knotIntervals.size(); ++i) {
				auto k1 = knotIntervals[i - 1];
				auto k2 = knotIntervals[i];
				
				auto spanH = find_span(host, k1);
				gather_weights(model, host, spanH, host_weights);

				auto positionsH = findGlobalPositions(host, spanH);
				ozp::quadrature::integrate(quad, ozp::quadrature::make_interval(k1, k2), [&](T u, T w){
					compute_nurbs_basis(host, spanH, u, host_weights, host_r);
					pointReverser(model, host, guest, spanH, u, spanG, uG);

					gather_weights(model, guest, spanG, guest_weights);
					compute_nurbs_basis(guest, spanG, uG, guest_weights, guest_r);
					mortar_l.noalias() = w * guest_r * host_r.transpose();
					auto positionsG = findGlobalPositions(guest, spanG);
					assemble_mortar_matrix(positionsH, positionsG, mortar_l, mortar);
				});
			}
		}

    template <template<class, unsigned int> class NurbsType, typename T, unsigned int DimG, typename PointReverserType>
		void mortar_integrator(const squanch::Model<T>& model, PointReverserType& pointReverser, const ozp::quadrature::Quadrature& quad, const NurbsType<T, 2>& host, const NurbsType<T, DimG>& guest, 
			Eigen::Matrix<T, -1, -1>& mortar)
    {
      const T small_number = 0.0000001;
			mortar.setZero();

      auto num_local_nodes_host = num_nodes_per_span(host.curves());
      auto num_local_nodes_guest = num_nodes_per_span(guest.curves());

      Eigen::Matrix<T, -1, -1> mortar_l(num_local_nodes_guest, num_local_nodes_host);
      Eigen::Matrix<T, -1, 1> host_weights(num_local_nodes_host);
      Eigen::Matrix<T, -1, 1> guest_weights(num_local_nodes_guest);
      Eigen::Matrix<T, -1, 1> host_r(num_local_nodes_host);
      Eigen::Matrix<T, -1, 1> guest_r(num_local_nodes_guest);

      NurbsType<T, DimG>::ParametricPointType uG;
      NurbsType<T, DimG>::SpanType spanG;
      NurbsType<T, 2>::SpanType spanH;
      NurbsType<T, 2>::ParametricPointType uH;

			auto&& knotsU = detail::findSmallestKnotIntervals(host.curves()[0].knots(), guest.curves()[0].knots());
			auto&& knotsV = detail::findSmallestKnotIntervals(host.curves()[1].knots(), guest.curves()[1].knots());

			for (auto i = 1u; i < knotsU.size(); ++i) {
				auto ku1 = knotsU[i - 1];
				auto ku2 = knotsU[i];
				for (auto j = 1u; j < knotsV.size(); ++j) {
					auto kv1 = knotsV[j - 1];
					auto kv2 = knotsV[j];
					spanH = find_span(host, std::array<T, 2>({ { ku1, kv1 } }));
					gather_weights(model, host, spanH, host_weights);
					auto positionsH = findGlobalPositions(host, spanH);

					ozp::quadrature::integrate(quad, ozp::quadrature::make_interval(ku1, kv1, ku2, kv2), [&](T ux, T uy, T w1, T w2){
						uH = { { ux, uy } };
						compute_nurbs_basis(host, spanH, uH, host_weights, host_r);
						pointReverser(model, host, guest, spanH, uH, spanG, uG);

						gather_weights(model, guest, spanG, guest_weights);
						compute_nurbs_basis(guest, spanG, uG, guest_weights, guest_r);

						mortar_l.noalias() = w1 * w2 * guest_r * host_r.transpose();
						auto positionsG = findGlobalPositions(guest, spanG);
						assemble_mortar_matrix(positionsH, positionsG, mortar_l, mortar);
					});
				}
			}
		}

		template <typename T, template<class, unsigned int> class NurbsType, unsigned int DimM, unsigned int DimS, typename MortarType> void create_connection_structure
			(const NurbsType<T, DimM>& master, const NurbsType<T, DimS>& slave, const MortarType& mortar, std::vector<std::pair<int, std::vector<std::pair<int, T>>>>& connection)
		{
			const T small_number = 1e-8;
      auto slave_ids = std::vector<int>(slave.cpi().begin(), slave.cpi().end());
      auto master_ids = std::vector<int>(master.cpi().begin(), master.cpi().end());

			connection.resize(mortar.rows());

			for (int k = 0; k < mortar.outerSize(); ++k) {
				for (MortarType::InnerIterator it(mortar, k); it; ++it) {
					auto coeff = it.value();
					if (coeff > small_number || coeff < -small_number) {
						auto&& conn = connection[it.row()];
						conn.first = slave_ids.data()[it.row()];
						conn.second.push_back(std::make_pair(master_ids.data()[it.col()], coeff));
					}
				}
			}
		}
	}

	template <typename T> void clear_out_small_numbers(const Eigen::Matrix<T, -1, -1>& mortar, T err, Eigen::SparseMatrix<T>& mortar_sparse)
	{
		mortar_sparse.resize(mortar.rows(), mortar.cols());
		std::vector<Eigen::Triplet<T>> trips;

		for (auto i = 0u; i < mortar.rows(); ++i) {
			for (auto j = 0u; j < mortar.cols(); ++j) {
				auto val = mortar(i,j);
				if (std::abs(val) > err) trips.emplace_back(i, j, val);
			}
		}
		mortar_sparse.setFromTriplets(trips.begin(), trips.end());
	}

	void solve_dense(Eigen::MatrixXd& k, Eigen::MatrixXd& f)
	{
		auto uplo = LAPACK_COL_MAJOR;
		auto n = (int)k.rows();
		auto nrhs = (int)f.cols();
		auto a = k.data();
		std::vector<int> ipiv(n);
		auto b = f.data();
		auto lda = (int)k.outerStride();
		auto ldb = (int)f.outerStride();
		int info = 0;

		dgesv(&n, &nrhs, a, &lda, ipiv.data(), b, &ldb, &info);

		if (info != 0) {
			throw std::runtime_error("Error in solving dense linear equation. MKL Error Code:" + std::to_string(info));
		}
	}

	/// <summary>
	/// Connects the slave patch to master patch using the mortar method
	/// </summary>
	/// <param name="model">The model.</param>
	/// <param name="master">The master.</param>
	/// <param name="slave">The slave.</param>
	/// <param name="quad">The quadrature to use.</param>
	/// <param name="mortar">The mortar.</param>
	/// <param name="slave_to_master">The slave_to_master connection data. [<slave_node_id, [<master_node_id, connection_coeff>]>]		
	/// </param>
	template <template<class, unsigned int> class NurbsType, typename T, unsigned int DimM, unsigned int DimS>
  void connect_with_mortar
    (const squanch::Model<T>& model, const NurbsType<T, DimM>& master, const NurbsType<T, DimS>& slave, const ozp::quadrature::Quadrature& quad, Eigen::Matrix<T, -1, -1>& mortar,
		std::vector<std::pair<int, std::vector<std::pair<int, T>>>>& slave_to_master)
  {
    const T small_number = 0.0000001;

    auto num_nodes_slave_g = slave.cpi().size();
    auto num_nodes_master_g = master.cpi().size();
		
		mortar.resize(num_nodes_slave_g, num_nodes_master_g);
		Eigen::Matrix<T, -1, -1> dsg(num_nodes_slave_g, num_nodes_slave_g);

    detail::IdentityPointReverser identityPR;
    detail::NewtonRaphsonPointReverser nrPR;
    detail::mortar_integrator(model, identityPR, quad, slave, slave, dsg);
    detail::mortar_integrator(model, nrPR, quad, master, slave, mortar);

		solve_dense(dsg, mortar);
		Eigen::SparseMatrix<double> mortar_sparse;
		clear_out_small_numbers(mortar, 1e-8, mortar_sparse);
		detail::create_connection_structure(master, slave, mortar_sparse, slave_to_master);
  }

	template <template<class, unsigned int> class NurbsType, typename T, unsigned int DimM, unsigned int DimS>
	void connect_with_mortar
		(const squanch::Model<T>& model, const NurbsType<T, DimM>& master, const NurbsType<T, DimS>& slave, const ozp::quadrature::Quadrature& quad,
		std::vector<std::pair<int, std::vector<std::pair<int, T>>>>& slave_to_master)
	{
		Eigen::Matrix<T, -1, -1> mortar;
		connect_with_mortar(model, master, slave, quad, mortar, slave_to_master);
	}

	template <template<class, unsigned int> class NurbsType, typename T, unsigned int DimM, unsigned int DimS>
	void connect_with_mortar_two_way
		(const squanch::Model<T>& model, const NurbsType<T, DimM>& master, const NurbsType<T, DimS>& slave, const ozp::quadrature::Quadrature& quad,
		std::vector<std::pair<int, std::vector<std::pair<int, T>>>>& slave_to_master, std::vector<std::pair<int, std::vector<std::pair<int, T>>>>& master_to_slave)
	{
		Eigen::Matrix<T, -1, -1> mortar;
		connect_with_mortar(model, master, slave, quad, mortar, slave_to_master);
		detail::create_connection_structure(slave, master, mortar.transpose(), master_to_slave);
	}
}