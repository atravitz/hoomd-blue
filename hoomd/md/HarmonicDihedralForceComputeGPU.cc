// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.


// Maintainer: dnlebard

/*! \file HarmonicDihedralForceComputeGPU.cc
    \brief Defines HarmonicDihedralForceComputeGPU
*/



#include "HarmonicDihedralForceComputeGPU.h"

namespace py = pybind11;
using namespace std;

/*! \param sysdef System to compute bond forces on
*/
HarmonicDihedralForceComputeGPU::HarmonicDihedralForceComputeGPU(std::shared_ptr<SystemDefinition> sysdef)
    : HarmonicDihedralForceCompute(sysdef)
    {
    // can't run on the GPU if there aren't any GPUs in the execution configuration
    if (!m_exec_conf->isCUDAEnabled())
        {
        m_exec_conf->msg->error() << "Creating a DihedralForceComputeGPU with no GPU in the execution configuration" << endl;
        throw std::runtime_error("Error initializing DihedralForceComputeGPU");
        }

    // allocate and zero device memory
    GPUArray<Scalar4> params(m_dihedral_data->getNTypes(),m_exec_conf);
    m_params.swap(params);

    m_tuner.reset(new Autotuner(32, 1024, 32, 5, 100000, "harmonic_dihedral", this->m_exec_conf));
    }

HarmonicDihedralForceComputeGPU::~HarmonicDihedralForceComputeGPU()
    {
    }

/*! \param type Type of the dihedral to set parameters for
    \param K Stiffness parameter for the force computation
    \param sign the sign of the cosine term
        \param multiplicity the multiplicity of the cosine term

    Sets parameters for the potential of a particular dihedral type and updates the
    parameters on the GPU.
*/
void HarmonicDihedralForceComputeGPU::setParams(unsigned int type, Scalar K, int sign, unsigned int multiplicity)
    {
    HarmonicDihedralForceCompute::setParams(type, K, sign, multiplicity);

    ArrayHandle<Scalar4> h_params(m_params, access_location::host, access_mode::readwrite);
    // update the local copy of the memory
    h_params.data[type] = make_scalar4(Scalar(K), Scalar(sign), Scalar(multiplicity), Scalar(0.0));
    }

/*! Internal method for computing the forces on the GPU.
    \post The force data on the GPU is written with the calculated forces

    \param timestep Current time step of the simulation

    Calls gpu_compute_harmonic_dihedral_forces to do the dirty work.
*/
void HarmonicDihedralForceComputeGPU::computeForces(unsigned int timestep)
    {
    // start the profile
    if (m_prof) m_prof->push(m_exec_conf, "Harmonic Dihedral");

    ArrayHandle<DihedralData::members_t> d_gpu_dihedral_list(m_dihedral_data->getGPUTable(), access_location::device,access_mode::read);
    ArrayHandle<unsigned int> d_n_dihedrals(m_dihedral_data->getNGroupsArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_dihedrals_ABCD(m_dihedral_data->getGPUPosTable(), access_location::device, access_mode::read);

    // the dihedral table is up to date: we are good to go. Call the kernel
    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
    BoxDim box = m_pdata->getGlobalBox();

    ArrayHandle<Scalar4> d_force(m_force,access_location::device,access_mode::overwrite);
    ArrayHandle<Scalar> d_virial(m_virial,access_location::device,access_mode::overwrite);
    ArrayHandle<Scalar4> d_params(m_params, access_location::device, access_mode::read);

    // run the kernel in parallel on all GPUs
    this->m_tuner->begin();
    gpu_compute_harmonic_dihedral_forces(d_force.data,
                                         d_virial.data,
                                         m_virial.getPitch(),
                                         m_pdata->getN(),
                                         d_pos.data,
                                         box,
                                         d_gpu_dihedral_list.data,
                                         d_dihedrals_ABCD.data,
                                         m_dihedral_data->getGPUTableIndexer().getW(),
                                         d_n_dihedrals.data,
                                         d_params.data,
                                         m_dihedral_data->getNTypes(),
                                         this->m_tuner->getParam());
    if(m_exec_conf->isCUDAErrorCheckingEnabled())
        CHECK_CUDA_ERROR();
    this->m_tuner->end();

    if (m_prof) m_prof->pop(m_exec_conf);
    }

void export_HarmonicDihedralForceComputeGPU(py::module& m)
    {
    py::class_<HarmonicDihedralForceComputeGPU, std::shared_ptr<HarmonicDihedralForceComputeGPU> >(m, "HarmonicDihedralForceComputeGPU", py::base<HarmonicDihedralForceCompute>())
    .def(py::init< std::shared_ptr<SystemDefinition> >())
    ;
    }
