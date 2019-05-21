/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "enums.hh"

std::ostream& operator<<(std::ostream& out, PrecondType preconditioner)
{
  switch (preconditioner) {
    case PrecondType::NONE: out << "NONE"; break;
    case PrecondType::JACOBI: out << "JACOBI"; break;
    case PrecondType::PARTIAL: out << "PARTIAL"; break;
    case PrecondType::ILUT: out << "ILUT"; break;
    case PrecondType::ILUK: out << "ILUK"; break;
    case PrecondType::ARMS: out << "ARMS"; break;
    case PrecondType::INVERSE: out << "INVERSE"; break;
    case PrecondType::ADDITIVESCHWARZ: out << "ADDITIVESCHWARZ"; break;
    case PrecondType::BOOMERAMG: out << "BOOMERAMG"; break;
    case PrecondType::EUCLID: out << "EUCLID"; break;
    case PrecondType::TRILINOSML: out << "TRILINOSML"; break;
    case PrecondType::SSOR: out << "SSOR"; break;
    case PrecondType::ICC0: out << "ICC0"; break;
    case PrecondType::ICC: out << "ICC"; break;
    case PrecondType::ILUKS: out << "ILUKS"; break;
    case PrecondType::HB: out << "HB"; break;
    case PrecondType::DIRECT: out << "DIRECT"; break;
    default: out << "(Unknown preconditioner key)"; break;
  }
  
  return out;
};


std::ostream& operator<<(std::ostream& out, IterateType iterative)
{
  switch (iterative) {
    case IterateType::CG: out << "CG"; break;
    case IterateType::BICGSTAB: out << "BICGSTAB"; break;
    case IterateType::GMRES: out << "GMRES"; break;
    case IterateType::PCG: out << "PCG"; break;
    case IterateType::APCG: out << "APCG"; break;
    case IterateType::SGS: out << "SGS"; break;
    default: out << "(Unknown iterative solver key)"; break;
  }
  
  return out;
};


std::ostream& operator<<(std::ostream& out, DirectType direct)
{
  switch (direct) {
    case DirectType::ANY:         out << "any one"; break;
    case DirectType::UMFPACK:     out << "UMFPACK"; break;
    case DirectType::PARDISO:     out << "PARDISO"; break;
    case DirectType::MUMPS:       out << "MUMPS"; break;
    case DirectType::SUPERLU:     out << "SUPERLU"; break;
    case DirectType::UMFPACK3264: out << "UMFPACK3264"; break;
    case DirectType::UMFPACK64:   out << "UMFPACK64"; break;
    default:                      out << "(Unknown direct solver key)"; break;
  }
  
  return out;
};

std::ostream& operator<<(std::ostream& out, MatrixProperties properties)
{
  switch (properties) {
    case MatrixProperties::GENERAL: out << "GENERAL"; break;
    case MatrixProperties::SYMMETRICSTRUCTURE: out << "SYMMETRICSTRUCTURE"; break;
    case MatrixProperties::SYMMETRIC: out << "SYMMETRIC"; break;
    case MatrixProperties::POSITIVEDEFINITE: out << "POSITIVEDEFINITE"; break;
    default: out << "(Unknown matrix properties key)"; break;
  }
  
  return out;
};

