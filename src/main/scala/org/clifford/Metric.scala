
package org.clifford

import scala.collection.mutable.Buffer

import breeze.linalg._, eigSym.EigSym

trait Metric {

  /**
    * @return true iff this metric is Euclidean
    * <p>There is no tolerance. The metric must be exactly Euclidean.
    */
  def isEuclidean(): Boolean

  /**
    * @return true iff this metric is anti-Euclidean
    * <p>There is no tolerance. The metric must be exactly Anti-Euclidean.
    */
  def isAntiEuclidean(): Boolean
  
  /**
    * @return true iff this metric non-zero entries has only only
    * on the diagonal of it's matrix. 
    *
    * <p>In other words: {e_i . e_j = 0} if (i != j)
    *
    * <p>There is no tolerance. Even the smallest non zero off-diagonal entry
    * will cause this method to return false.
    */
  def isDiagonal(): Boolean
    
  /**
    * Returns the symmetric, real matrix representing the metric.
    *
    * <p>Each entry M[i, j] holds the value of e_i . e_j
    *
    * <p>This matrix may or may not be created from scratch on each
    * call (so don't count on it being efficient, although it is easy
    * to 'cache' the result of the first call to this method when efficiency
    * is required).
     *
    * <p>Do not modify the returned matrix.
    *
    */
  def getInnerProductMatrix(): DenseMatrix[Double]
  //def getInnerProductMatrix(): cern.colt.matrix.DoubleMatrix2D
    
  /**
    * @return the value of the inner product: e_idx .e_idx
    */
  def getBasisVectorIp(idx: Int): Double
    
  /**
    * @return the value of the inner product: e_idx1 .e_idx2
    */
  def getBasisVectorIp(idx1: Int, idx2: Int): Double

  /** transforms a to the metric basis (A must be on eigenbasis) */
  def toMetricBasis(a: BasisBlade): Seq[BasisBlade]

  /** transforms A to the metric basis (a must be on eigenbasis) */
  def toMetricBasis(a: Seq[BasisBlade]): Seq[BasisBlade]

  /** transforms a to the eigen basis (a must be on metric basis) */
  def toEigenbasis(a: BasisBlade): Seq[BasisBlade]

  /** transforms A to the eigen basis (A must be on metric basis) */
  def toEigenbasis(a: Seq[BasisBlade]): Seq[BasisBlade]
    
  /**
    * Returns the eigenvalue decomposition of the symmetric, real matrix representing the metric.
    *
    * <p>This EigenvalueDecomposition may or may not be created from scratch on each
    * call (so don't count on it being efficient, although it is easy
    * to 'cache' the result of the first call to this method when efficiency
    * is required).
    *
    * <p>Do not modify the returned EigenvalueDecomposition.
    */
  //def getInnerProductMatrixEigenvalueDecomposition(): cern.colt.matrix.linalg.EigenvalueDecomposition
  def getInnerProductMatrixEigenvalueDecomposition():
      (DenseVector[Double], DenseMatrix[Double])

  val eigenMetric: DenseVector[Double]
}

sealed trait InnerProductType
case object LeftContraction extends InnerProductType
case object RightContraction extends InnerProductType
case object HestenesInnerProduct extends InnerProductType
case object ModifiedHestenesInnerProduct extends InnerProductType

object Metric {

  def apply(matrix: Array[Array[Double]]): Metric = 
    apply(DenseMatrix(matrix: _*))

  def apply(matrix: DenseMatrix[Double]): Metric = new MetricImpl(matrix)
}

/** 
  * @arg matrix the metric (symmetric matrix) 
  */
class MetricImpl(val matrix: DenseMatrix[Double]) extends Metric {

  /** the eigenvectors matrix & eigenvalues of m_matrix */
  private val eigen = eigSym(matrix)
  val eigenDecom: DenseMatrix[Double] = eigen.eigenvectors
  val eigenMetric: DenseVector[Double] = eigen.eigenvalues

  /** inverse of the eigenmatrix */
  //protected DoubleMatrix2D m_invEigMatrix;
  val invEigen: DenseMatrix[Double] = eigenDecom.t

  def isDiagonal(): Boolean =
    matrix.forall { (loc, v) => (loc._1 == loc._2 || v == 0.0) }
  def isEuclidean(): Boolean = isDiagonal() && diag(matrix).forall { _ == 1.0 }
  def isAntiEuclidean(): Boolean =
    isDiagonal() && diag(matrix).forall { _ == -1.0 }

  def getBasisVectorIp(idx: Int): Double = getBasisVectorIp(idx, idx)

  def getBasisVectorIp(idx1: Int, idx2: Int): Double = matrix(idx1, idx2)

  def getInnerProductMatrix(): DenseMatrix[Double] = matrix

  def getInnerProductMatrixEigenvalueDecomposition(): 
      (DenseVector[Double], DenseMatrix[Double]) = (eigenMetric, eigenDecom)

  /** transforms a to the eigen basis (a must be on metric basis) */
  def toEigenbasis(a: BasisBlade): Seq[BasisBlade] = transform(a, invEigen)

  /** transforms A to the eigen basis (A must be on metric basis) */
  def toEigenbasis(a: Seq[BasisBlade]): Seq[BasisBlade] =
    BasisBlade.simplify(a.flatMap { toEigenbasis(_) })

  /** transforms a to the metric basis (A must be on eigenbasis) */
  def toMetricBasis(a: BasisBlade): Seq[BasisBlade] = transform(a, eigenDecom)

  /** transforms A to the metric basis (a must be on eigenbasis) */
  def toMetricBasis(a: Seq[BasisBlade]): Seq[BasisBlade] =
    BasisBlade.simplify(a.flatMap { toMetricBasis(_) })

  protected def transform(
    a: BasisBlade, m: DenseMatrix[Double]): Seq[BasisBlade] = {

    var ab = Seq[BasisBlade](new BasisBlade(a.scale))
    var i = 0
    var b = a.bitmap

    while (b != 0) {
      if ((b & 1) != 0) {
        val newAb = Buffer[BasisBlade]()
        for (j <- (0 until m.rows)) {
          if (m(j, i) != 0.0) {
            val diagValue = m(j, i)
            newAb ++= ab.map { blade =>
              blade ^ (new BasisBlade((1 << j), diagValue)) }
          }
        }
        ab = newAb.toSeq
      }

      b = b >> 1
      i = i + 1
    }

    ab
  }
}
