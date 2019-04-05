
package org.clifford

import org.scalatest.{ FlatSpec, Matchers }

class BasisBladeSpace extends FlatSpec with Matchers with CustomMatchers {

  behavior of "BasisBlade"

  it should "rounding operations should work" in {

    import BasisBlade._

    implicit val tol = 0.000000000001
    val a = Array(1.0001, 1.49, 2.02, -0.51, -10000.1)

    roundDouble(a(0), 0.07, 0.07) should beApproximately(0.98)
    roundDouble(a(1), 0.07, 0.07) should beApproximately(1.47)
    roundDouble(a(2), 0.07, 0.07) should beApproximately(2.03)
    roundDouble(a(3), 0.07, 0.07) should beApproximately(-0.49)
    roundDouble(a(4), 0.07, 0.07) should beApproximately(-10000.13)
  }

  it should "parse implicit metrics" in {

    val m = Array(
      Array(1.0, 0.0, 0.0), 
      Array(0.0, 0.0, -1.0), 
      Array(0.0, -1.0, 1.0))

    implicit val metric: Option[Metric] = Some(m)

    // TODO more tests here
  }

  it should "complete basic operations" in {
  
    val e1 = new BasisBlade(1)
    val no = new BasisBlade(2)
    val ni = new BasisBlade(4)

    val m: Array[Array[Double]] = Array(
      Array(1.0, 0.0, 0.0, 0.0, 0.0),
      Array(0.0, 1.0, 0.0, 0.0, 0.0),
      Array(0.0, 0.0, 1.0, 0.0, 0.0),
      Array(0.0, 0.0, 0.0, 0.0, -1.0),
      Array(0.0, 0.0, 0.0, -1.0, 0.0)
    )

    //val m = new Metric(m)

    val t = System.currentTimeMillis()
    for (i <- 0 until 32; j <- 0 until 32) {
      val a = new BasisBlade(i)
      val b = new BasisBlade(j)
      //val r = BasisBlade.gp(a, b, m).round(1.0, 0.0001)
      val r = a.innerProduct(b).round(1.0, 0.0001)
      //(a * b).round(1.0, 0.0001) should be (new BasisBlade(1.0))
      //println(s"rounded $a * $b = $r")

      val r1 = a * b //BasisBlade.gp(a, b, m)
      //r1.grade() should be(2)
      //println(s"$a * $b = $r1")
    }
       
    println("time: " + (System.currentTimeMillis() - t))
  }
}

