//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------
#include "MantidPythonInterface/api/WorkspaceToNumpy.h"
#include "MantidPythonInterface/kernel/NumpyConverters.h"
#include "MantidAPI/MatrixWorkspace.h"

#include <boost/python/extract.hpp>

// See http://docs.scipy.org/doc/numpy/reference/c-api.array.html#PY_ARRAY_UNIQUE_SYMBOL
#define PY_ARRAY_UNIQUE_SYMBOL KERNEL_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/ndarrayobject.h>

namespace Mantid
{
  namespace PythonInterface
  {
    namespace Numpy
    {
      using Mantid::API::MatrixWorkspace_sptr;
      using Mantid::API::MatrixWorkspace;
      namespace bpl = boost::python;

      // ----------------------------------------------------------------------------------------------------------
      namespace
      {
        /// Which data field are we extracting
        enum DataField { XValues = 0, YValues = 1, EValues= 2 };

        /**
         * Helper method for extraction to numpy.
         * @param workspace :: A pointer to the workspace that contains the data
         * @param field :: Which field should be extracted
         * @param start :: The index in the workspace to start at when reading the data
         * @param endp1 :: One past the end index in the workspace to finish at when reading the data (similar to .end() for STL)
         *
         */
        PyArrayObject *cloneArray(MatrixWorkspace_sptr workspace,
                                   DataField field, const size_t start, const size_t endp1)
        {
          const size_t numHist = endp1 - start;
          npy_intp stride(0);

          // Find out which function we need to call to access the data
          typedef const MantidVec & (MatrixWorkspace::*ArrayAccessFn)(const size_t) const;
          ArrayAccessFn dataAccesor;
          if( field == XValues )
          {
            stride = workspace->readX(0).size();
            dataAccesor = &MatrixWorkspace::readX;
          }
          else
          {
            stride = workspace->blocksize();
            if(field == YValues) dataAccesor = &MatrixWorkspace::readY;
            else dataAccesor = &MatrixWorkspace::readE;
          }
          npy_intp arrayDims[2] = {static_cast<npy_intp>(numHist), stride};
          PyArrayObject *nparray = (PyArrayObject *)PyArray_NewFromDescr(&PyArray_Type,
                                                                         PyArray_DescrFromType(NPY_DOUBLE),
                                                                         2, // rank 2
                                                                         arrayDims, // Length in each dimension
                                                                         NULL, NULL,
                                                                         0, NULL);
          double *dest = (double*)PyArray_DATA(nparray); // HEAD of the contiguous numpy data array
          for(size_t i = start; i < endp1; ++i)
          {
            const MantidVec & src = ((*workspace).*(dataAccesor))(i);
            std::copy(src.begin(), src.end(), dest);
            dest += stride; // Move the ptr to the start of the next 1D array
          }
          return nparray;
        }
      }


      // -------------------------------------- Read-only arrays---------------------------------------------------
      /*
       * Create a numpy wrapper around the original X values at the given index
       * @param self :: A pointer to a PyObject representing the calling object
       * @param index :: The index into the workspace
       */
      PyObject *wrapX(PyObject *self, const size_t index)
      {
        MatrixWorkspace_sptr workspace = bpl::extract<MatrixWorkspace_sptr>(self);
        PyArrayObject *nparray = (PyArrayObject*)wrapWithReadOnlyNumpy(workspace->readX(index));
        return (PyObject*)nparray;
      }
      /*
       * Create a numpy wrapper around the original Y values at the given index
       * @param self :: A pointer to a PyObject representing the calling object
       * @param index :: The index into the workspace
       */
      PyObject *wrapY(PyObject *self, const size_t index)
      {
        MatrixWorkspace_sptr workspace = bpl::extract<MatrixWorkspace_sptr>(self);
        PyArrayObject *nparray = (PyArrayObject*)wrapWithReadOnlyNumpy(workspace->readY(index));
        return (PyObject*)nparray;
      }

      /*
       * Create a numpy wrapper around the original E values at the given index
       * @param self :: A pointer to a PyObject representing the calling object
       * @param index :: The index into the workspace
       */
      PyObject *wrapE(PyObject *self, const size_t index)
      {
        MatrixWorkspace_sptr workspace = bpl::extract<MatrixWorkspace_sptr>(self);
        PyArrayObject *nparray = (PyArrayObject*)wrapWithReadOnlyNumpy(workspace->readE(index));
        return (PyObject*)nparray;
      }


      // -------------------------------------- Cloned arrays---------------------------------------------------
      /* Create a numpy array from the X values of the given workspace reference
       * This acts like a python method on a Matrixworkspace object
       * @param self :: A pointer to a PyObject representing the calling object
       * @return A 2D numpy array created from the X values
       */
      PyObject * cloneX(PyObject * self)
      {
        MatrixWorkspace_sptr workspace = bpl::extract<MatrixWorkspace_sptr>(self);
        return (PyObject*)cloneArray(workspace, XValues, 0, workspace->getNumberHistograms());
      }
      /* Create a numpy array from the Y values of the given workspace reference
       * This acts like a python method on a Matrixworkspace object
       * @param self :: A pointer to a PyObject representing the calling object
       * @return A 2D numpy array created from the Y values
       */
      PyObject * cloneY(PyObject * self)
      {
        MatrixWorkspace_sptr workspace = bpl::extract<MatrixWorkspace_sptr>(self);
        return (PyObject*)cloneArray(workspace, YValues, 0, workspace->getNumberHistograms());
      }

      /* Create a numpy array from the E values of the given workspace reference
       * This acts like a python method on a Matrixworkspace object
       * @param self :: A pointer to a PyObject representing the calling object
       * @return A 2D numpy array created from the E values
       */
      PyObject * cloneE(PyObject * self)
      {
        MatrixWorkspace_sptr workspace = bpl::extract<MatrixWorkspace_sptr>(self);
        return (PyObject*)cloneArray(workspace, EValues, 0, workspace->getNumberHistograms());
      }

    }
  }
}



