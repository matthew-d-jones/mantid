#ifndef MANTID_MDEVENTS_COORDTRANSFORM_H_
#define MANTID_MDEVENTS_COORDTRANSFORM_H_
    
#include "MantidKernel/System.h"
#include "MantidGeometry/MDGeometry/IMDDimension.h"
#include "MantidGeometry/Math/Matrix.h"


namespace Mantid
{
namespace MDEvents
{
  /** Template for shortening template coordinate transform classes. */
  #define TCT template<size_t inD, size_t outD>

  /** Generic class to transform from M input dimensionsto N output dimensions.
   *
   * The types of conversions to account for are:
   * * Simple rotation matrix
   * * Affine Transformation = linear transform such as a rotation + a translation
   * * Projection into lower dimensions, for example taking a 2D slice out of 3D data.
   * 
   * This class could be subclassed in order to handle non-linear transforms (though
   * making the apply() method virtual would disallow inlining = slowdown).
   *
   * @author Janik Zikovsky
   * @date 2011-04-14 10:03:55.944809
   */
  class DLLExport CoordTransform 
  {
  public:
    CoordTransform(const size_t inD, const size_t outD);
    ~CoordTransform();

    void addTranslation(const CoordType * translationVector);
    Mantid::Geometry::Matrix<CoordType> getMatrix() const;
    void setMatrix(const Mantid::Geometry::Matrix<CoordType> newMatrix);

    //----------------------------------------------------------------------------------------------
    /** Apply the coordinate transformation
     *
     * @param inputVector :: fixed-size array of input coordinates, of size inD
     * @param outVector :: fixed-size array of output coordinates, of size outD
     */
    void apply(const CoordType * inputVector, CoordType * outVector)
    {
      // For each output dimension
      for (size_t out = 0; out < outD; ++out)
      {
        //Cache the row pointer to make the matrix access a bit faster
        CoordType * rawMatrixRow = rawMatrix[out];
        CoordType outVal = 0.0;
        size_t in;
        for (in = 0; in < inD; ++in)
          outVal += rawMatrixRow[in] * inputVector[in];

        // The last input coordinate is "1" always (made homogenous coordinate out of the input x,y,etc.)
        outVal += rawMatrixRow[in];
        // Save in the output
        outVector[out] = outVal;
      }
    }



  protected:
    /// Input number of dimensions
    size_t inD;

    /// Output number of dimensions
    size_t outD;

    /** Affine Matrix to perform the transformation. The matrix has inD+1 columns, outD+1 rows.
     * By using an affine, translations and rotations (or other linear transforms) can be
     * combined by simply multiplying the matrices.
     */
    Mantid::Geometry::Matrix<CoordType> affineMatrix;

    /// Raw pointer to the same underlying matrix as affineMatrix.
    CoordType ** rawMatrix;

    void copyRawMatrix();
  };


} // namespace Mantid
} // namespace MDEvents

#endif  /* MANTID_MDEVENTS_COORDTRANSFORM_H_ */
