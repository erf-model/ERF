#ifndef ERF_TIME_INTERP_DATA_H_
#define ERF_TIME_INTERP_DATA_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

struct TimeInterpolatedData {
    TimeInterpolatedData() {}

    TimeInterpolatedData(amrex::Real time) : m_time(time) {}

    ~TimeInterpolatedData() {
        clear();
    }

    void set_time(amrex::Real time) {
        m_time = time;
    }

    amrex::Real get_time() { return m_time; }

    void add_var(amrex::MultiFab* var_data, int own_data) {
        m_data.push_back(var_data);
        m_owns_data.push_back(own_data);
    }

    amrex::MultiFab& get_var(int var_idx) {
        return *m_data[var_idx];
    }

    int num_vars() { return m_data.size(); }

    void clear() {
        // clear all data owned by this object
        for (int i = 0; i < num_vars(); ++i) {
            if (m_owns_data[i] == owning) {
                m_data[i]->clear();
                delete m_data[i];
            }
        }

        // clear internal memory
        m_data.clear();
        m_owns_data.clear();
        m_time = 0.;
    }

    const int owning = 1;
    const int non_owning = 0;

private:
    amrex::Vector<amrex::MultiFab*> m_data;
    amrex::Vector<int> m_owns_data;
    amrex::Real m_time;
};

#endif
