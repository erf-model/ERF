#ifndef MATERIALPROPERTIES_H
#define MATERIALPROPERTIES_H

#include <string>
#include <AMReX_REAL.H>

/*! Base class for material properties */
class MaterialProperties {
    public:

        MaterialProperties () { }

        virtual ~MaterialProperties () = default;

        virtual amrex::Real density() const = 0;

        virtual const std::string& name() const = 0;

    protected:

    private:
};

#endif
