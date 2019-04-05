
package org.clifford

import org.scalatest.matchers.{ MatchResult, Matcher }

trait CustomMatchers {

  class DoubleApproximatelyEqualTo(expectedValue: Double)(
    implicit tolerance: Double) extends Matcher[Double] {
    def apply(left: Double) = {
      MatchResult(
        (expectedValue - tolerance) <= left & left <= (expectedValue + tolerance),
        s"$left was not approximately equal to $expectedValue ($tolerance)",
        s"$left was approximately equal to $expectedValue ($tolerance)")
    }
  }

  def beApproximately(expectedValue: Double)(implicit tolerance: Double) =
    new DoubleApproximatelyEqualTo(expectedValue)
}
