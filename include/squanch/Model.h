#pragma once
#include <unordered_map>
#include <memory>
#include <squanch/PointContainerTypes.h>
#include <squanch/PatchBase.h>
#include <squanch/Nurbs.h>

namespace squanch {

    namespace detail {
        template <typename T> bool is_same(typename const squanch::PointContainerTypes<T>::point_type& p1, typename const squanch::PointContainerTypes<T>::point_type& p2, T epsilon)
        {
            return (p1 - p2).squaredNorm() < epsilon;
        }

/*        template <typename T, typename FunType> void for_each_point(const squanch::Model<T>& model, const std::vector<squanch::PatchBase*>& patches, FunType fun)
        {
            for (auto&& pb : patches) {
                auto&& cpi = pb->get_node_ids();
                for (auto i = 0u; i < cpi.size(); ++i) {
                    auto id = cpi.data()[i];
                    auto&& p = model.get_point(id);
                    fun(id, p);
                }
            }
        }*/

        //inline void replace_cpi(std::vector<squanch::PatchBase*>& patches, const int find_id, const int replace_id)
        //{
        //    for (auto&& pb : patches) {
        //        auto&& cpi = pb->get_node_ids();
        //        for (auto i = 0u; i < cpi.size(); ++i) {
        //            auto& id = cpi.data()[i];
        //            if (find_id == id) id = replace_id;
        //        }
        //    }
        //}

    }

    template <typename T> class Model {
    public:
        Model() : m_last_patch_id(1) {}

        const typename std::unordered_map<int, std::unique_ptr<squanch::PatchBase>>& patches() const { return m_patches; }

        template <int NDIM> squanch::Nurbs<T, NDIM>& new_nurbs()
        {
            auto p = new squanch::Nurbs<T, NDIM>(m_last_patch_id++);
            std::unique_ptr<squanch::PatchBase> nurbsp(p);
            m_patches[nurbsp->id()] = std::move(nurbsp);
            return *p;
        }

        void delete_patch(squanch::PatchBase* pb) { m_patches.erase(pb->id()); }
        squanch::PatchBase* get_patch(int id)
        {
            auto iter = m_patches.find(id);
            if (iter != m_patches.end()) return iter->second.get();
            else return nullptr;
        }

        std::pair<int, typename squanch::PointContainerTypes<T>::point_type&> new_point()
        {
            auto id = m_control_points.size();
            squanch::PointContainerTypes<T>::point_type p;
            m_control_points.push_back(p);
            return std::make_pair(id, std::ref(m_control_points.back()));
        }

        std::pair<int, typename squanch::PointContainerTypes<T>::point_type&>  new_point(T x, T y, T z, T w)
        {
            auto id = m_control_points.size();
            squanch::PointContainerTypes<T>::point_type p(x, y, z, w);
            m_control_points.push_back(p);
            return std::make_pair(id, std::ref(m_control_points.back()));
        }

        typename squanch::PointContainerTypes<T>::point_type& get_point(int id)
        {
            return m_control_points[id];
        }

        typename const squanch::PointContainerTypes<T>::point_type& get_point(int id) const
        {
            return m_control_points[id];
        }

        /*	void clear_unused_points(const std::vector<int>& used_point_ids)
            {
            squanch::PointContainerTypes<T>::map_type used_points;
            for (auto id : used_point_ids) {
            used_points[id] = m_control_points[id];
            }

            m_control_points = used_points;
            }
            */
        void clear()
        {
            m_control_points.clear();
            m_patches.clear();
            m_last_patch_id = 1;
        }

        size_t num_nodes() const { return m_control_points.size(); }

    private:
        typename squanch::PointContainerTypes<T>::vector_type m_control_points;
        std::unordered_map<int, std::unique_ptr<squanch::PatchBase>> m_patches;
        int m_last_patch_id;
    };
}