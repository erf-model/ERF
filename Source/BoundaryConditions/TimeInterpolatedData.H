#ifndef ERF_TIME_INTERP_DATA_H_
#define ERF_TIME_INTERP_DATA_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

struct TimeInterpolatedData {
    TimeInterpolatedData() = default;

    TimeInterpolatedData(amrex::Real time) : m_time(time) {}

    ~TimeInterpolatedData() {
        clear();
    }

    // Declare a default move constructor so we ensure the destructor is
    // not called when we return an object of this class by value
    TimeInterpolatedData(TimeInterpolatedData&&)  noexcept = default;

    // Declare a default move assignment operator
    TimeInterpolatedData& operator=(TimeInterpolatedData&& other)  noexcept = default;

    // Delete the copy constructor and copy assignment operators because
    // the integrator allocates internal memory that is best initialized
    // from scratch when needed instead of making a copy.

    // Delete the copy constructor
    TimeInterpolatedData(const TimeInterpolatedData& other) = delete;
    //
    // Delete the copy assignment operator
    TimeInterpolatedData& operator=(const TimeInterpolatedData& other) = delete;

    void set_time(amrex::Real time) {
        m_time = time;
    }

    [[nodiscard]] amrex::Real get_time() const { return m_time; }

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
            if (m_owns_data[i] == 1) {
                m_data[i]->clear();
                delete m_data[i];
            }
        }

        // clear internal memory
        m_data.clear();
        m_owns_data.clear();
        m_time = 0.;
    }

private:
    amrex::Vector<amrex::MultiFab*> m_data;
    amrex::Vector<int> m_owns_data;
    amrex::Real m_time;
};

#endif
