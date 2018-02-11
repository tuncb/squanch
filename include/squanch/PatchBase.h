#pragma once

namespace squanch {

  enum class PatchType {
    nurbs1d, nurbs2d, nurbs3d
  };

  class PatchBase {
  public:
    PatchBase(int id, PatchType pt) : m_type(pt), m_id(id) {}

    inline int         id() const { return m_id; }
    inline PatchType type() const { return m_type; }

  private:
    int m_id;
    PatchType m_type;
  };

}