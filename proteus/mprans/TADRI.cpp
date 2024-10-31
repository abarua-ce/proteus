#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "TADRI.h"

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
using proteus::TADRI_base;

PYBIND11_MODULE(cTADRI, m)
{
    xt::import_numpy();

    py::class_<TADRI_base>(m, "cTADRI_base")
        .def(py::init(&proteus::newTADRI))
        .def("calculateResidual", &TADRI_base::calculateResidual)
        .def("calculateJacobian", &TADRI_base::calculateJacobian)
        .def("FCTStep", &TADRI_base::FCTStep);
}
