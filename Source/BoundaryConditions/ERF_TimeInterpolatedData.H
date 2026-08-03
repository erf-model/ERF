#ifndef ERF_TIME_INTERP_DATA_H_
#define ERF_TIME_INTERP_DATA_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

struct TimeInterpolatedData {
    /**
     * Construct empty time-interpolated data storage.
     */
    TimeInterpolatedData() = default;

    /**
     * Construct storage associated with a time.
     *
     * @param[in] time time associated with the data
     */
    TimeInterpolatedData(double time) : m_time(time) {}

    /**
     * Destroy owned data and clear internal storage.
     */
    ~TimeInterpolatedData() {
        clear();
    }

    /**
     * Default move constructor.
     */
    TimeInterpolatedData(TimeInterpolatedData&&)  noexcept = default;

    /**
     * Default move-assignment operator.
     *
     * @param[in] other source object
     */
    TimeInterpolatedData& operator=(TimeInterpolatedData&& other)  noexcept = default;

    /**
     * Deleted copy constructor.
     *
     * @param[in] other source object
     */
    TimeInterpolatedData(const TimeInterpolatedData& other) = delete;

    /**
     * Deleted copy-assignment operator.
     *
     * @param[in] other source object
     */
    TimeInterpolatedData& operator=(const TimeInterpolatedData& other) = delete;

    /**
     * Set the time associated with the data.
     *
     * @param[in] time time associated with the data
     */
    void set_time(double time) {
        m_time = time;
    }

    /**
     * Return the time associated with the data.
     */
    [[nodiscard]] double get_time() const { return m_time; }

    /**
     * Add a MultiFab pointer to the managed list.
     *
     * @param[in] var_data MultiFab pointer to store
     * @param[in] own_data whether this object owns and deletes var_data
     */
    void add_var(amrex::MultiFab* var_data, int own_data) {
        m_data.push_back(var_data);
        m_owns_data.push_back(own_data);
    }

    /**
     * Return one stored MultiFab.
     *
     * @param[in] var_idx index of the stored MultiFab
     */
    amrex::MultiFab& get_var(int var_idx) {
        return *m_data[var_idx];
    }

    /**
     * Return the number of stored MultiFabs.
     */
    int num_vars() { return m_data.size(); }

    /**
     * Clear stored data, deleting entries owned by this object.
     */
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
        m_time = 0.0;
    }

private:
    amrex::Vector<amrex::MultiFab*> m_data;
    amrex::Vector<int> m_owns_data;
    double m_time;
};

#endif
