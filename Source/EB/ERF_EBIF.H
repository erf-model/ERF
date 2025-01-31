#ifndef ERF_EB_IF_H_
#define ERF_EB_IF_H_

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Union of a list (std::vector) of the same kind of implicit function.         *
 *                                                                              *
 ********************************************************************************/
template <class F>
class UnionListIF {

public:

    UnionListIF (const Vector<F> & a_ifs)
        : m_ifs(a_ifs)
    {
        m_empty = ! a_ifs.empty();
    }

    ~UnionListIF () {}

    UnionListIF (const UnionListIF & rhs) = default;
    UnionListIF (UnionListIF && rhs)      = default;
    UnionListIF & operator= (const UnionListIF & rhs) = default;
    UnionListIF & operator= (UnionListIF && rhs)      = default;

    bool is_empty () const {return m_empty;}

    Real operator() (const RealArray & p) const {

        // NOTE: this assumes that m_ifs is not empty
        Real vmax=m_ifs[0](p);
        for (int i=1; i< m_ifs.size(); i++) {
            Real vcur = m_ifs[i](p);
            if ( vmax < vcur ) vmax = vcur;
        }

        return vmax;
    }

private:

    Vector<F> m_ifs;
    bool m_empty;
};

#endif
