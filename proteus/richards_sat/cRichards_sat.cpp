#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "Richards_sat.h"

#if defined(__GNUC__) && !defined(__clang__)
    namespace workaround
    {
        inline void define_allocators()
        {
            std::allocator<int> a0;
            std::allocator<double> a1;
        }
    }
#endif

namespace py = pybind11;
using proteus::Richards_sat_base;

PYBIND11_MODULE(cRichards_sat, m)
{
    xt::import_numpy();

    py::class_<Richards_sat_base>(m, "cRichards_sat_base")
        .def(py::init(&proteus::newRichards_sat))
        .def("calculateResidual", &Richards_sat_base::calculateResidual)
        .def("calculateJacobian", &Richards_sat_base::calculateJacobian)
        .def("invert", &Richards_sat_base::invert)
        .def("FCTStep", &Richards_sat_base::FCTStep)
        .def("kth_FCT_step", &Richards_sat_base::kth_FCT_step)
        .def("calculateResidual_entropy_viscosity", &Richards_sat_base::calculateResidual_entropy_viscosity)
        .def("calculateMassMatrix", &Richards_sat_base::calculateMassMatrix);
}
