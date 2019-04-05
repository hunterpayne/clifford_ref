
package org.clifford

object MathU {

  /**
    * @return the hyperbolic cosine of @param x
    *
    * Computed as 0.5 (e^x + e^-x)
    */
  def cosh(x: Double): Double = 0.5 * (Math.exp(x) + Math.exp(-x))

  /**
    * @return the hyperbolic cosine of @param x
    *
    * Computed as 0.5 (e^x - e^-x)
    */
  def sinh(x: Double): Double = 0.5 * (Math.exp(x) - Math.exp(-x))
}
