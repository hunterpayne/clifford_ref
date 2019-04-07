
package org.clifford

import scala.collection.mutable.Buffer

object BasisBlade {

  /**
    * @return the number of 1 bits in <code>i</code>
    */
  def bitCount(arg: Int): Int = {
    var i = arg
    // Note that unsigned shifting (>>>) is not required.
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    i = (i + (i >> 4)) & 0x0F0F0F0F;
    i = i + (i >> 8);
    i = i + (i >> 16);
    i & 0x0000003F
  }

  def canonicalReorderingSign(a: Int, b: Int): Double = {
    // Count the number of basis vector flips required to
    // get a and b into canonical order.
    var a1: Int = a >> 1
    var sum: Int = 0

    while (a1 != 0) {
      sum = sum + bitCount(a1 & b)
      a1 >>= 1;
    }

    // even number of flips -> return 1
    // odd number of flips -> return -1
    if ((sum & 1) == 0) 1.0 else -1.0
  }

  def outerProduct(a: BasisBlade, b: BasisBlade): BasisBlade =
    // if outer product: check for independence
    if ((a.bitmap & b.bitmap) != 0) {
      new BasisBlade(0.0)
    } else {
      // compute the bitmap
      val bitmap: Int = a.bitmap ^ b.bitmap

      // compute the sign change due to reordering
      val sign: Double = canonicalReorderingSign(a.bitmap, b.bitmap)

      // return result
      BasisBlade(bitmap, sign * a.scale * b.scale)
    }

  /**
    * Computes the geometric product of two basis blades
    * @param m is an array of doubles giving the metric for each basis vector.
    */
  def geometricProduct(a: BasisBlade, b: BasisBlade): BasisBlade = {
    // compute the bitmap
    val bitmap: Int = a.bitmap ^ b.bitmap

    // compute the sign change due to reordering
    val sign: Double = canonicalReorderingSign(a.bitmap, b.bitmap)

    // return result
    BasisBlade(bitmap, sign * a.scale * b.scale)
  }

  def geometricProduct(
    a: BasisBlade, b: BasisBlade, arr: Array[Double]): BasisBlade = {
    // compute the geometric product in Euclidean metric:
    val result = geometricProduct(a, b)

    // compute the meet (bitmap of annihilated vectors):
    var bitmap = a.bitmap & b.bitmap

    var newScale = result.scale

    // change the scale according to the metric
    var i: Int = 0
    while (bitmap != 0) {
      if ((bitmap & 1) != 0) newScale = newScale * arr(i)
      i = i + 1
      bitmap = bitmap >> 1
    }

    BasisBlade(result.bitmap, newScale)
  }


  /**
    * Computes the geometric product of two basis blades in limited 
    * non-Euclidean metric.
    * @param m is an array of doubles giving the metric for each basis vector.
    */
  def geometricProductM(a: BasisBlade, b: BasisBlade)(
    implicit metric: Metric): Seq[BasisBlade] = {
    val aB = metric.toEigenbasis(a)
    val bB = metric.toEigenbasis(b)
    val metricArr = metric.eigenMetric.data

    metric.toMetricBasis(simplify(
      for { c <- aB; d <- bB } yield geometricProduct(c, d, metricArr)
    ))
  }

  /**
    * @return the geometric product of two basis blades
    * @param tpe gives the type of inner product:
    * LEFT_CONTRACTION,RIGHT_CONTRACTION, HESTENES_INNER_PRODUCT or 
    * MODIFIED_HESTENES_INNER_PRODUCT.
    */
  def innerProduct(a: BasisBlade, b: BasisBlade, tpe: InnerProductType)(
    implicit m: Option[Metric] = None) =
    innerProductFilterSeq(a.grade(), b.grade(), (a * b).blades.toSeq, tpe).
      headOption.getOrElse(new BasisBlade())

  private def innerProductFilterSeq(
    ga: Int, gb: Int, r: Seq[BasisBlade], 
    tpe: InnerProductType)(implicit m: Option[Metric] = None): Seq[BasisBlade] =
    r.map { r1 => innerProductFilter(ga, gb, r1, tpe) }.filter {
      _.scale != 0.0 }

  // *!*HTML_TAG*!* inner_product_filter
  /**
    * Applies the rules to turn a geometric product into an inner product
    * @param ga Grade of argument 'a'
    * @param gb Grade of argument 'b'
    * @param r the basis blade to be filter
    * @param type the type of inner product required:
    * LEFT_CONTRACTION,RIGHT_CONTRACTION, HESTENES_INNER_PRODUCT or 
    * MODIFIED_HESTENES_INNER_PRODUCT
    * @return Either a 0 basis blade, or 'r'
    */
  private def innerProductFilter(
    ga: Int, gb: Int, r: BasisBlade, tpe: InnerProductType): BasisBlade =
    tpe match {
      case LeftContraction =>
        if ((ga > gb) || (r.grade() != (gb - ga))) new BasisBlade()
        else r
      case RightContraction =>
        if ((ga < gb) || (r.grade() != (ga-gb))) new BasisBlade()
        else r
      case HestenesInnerProduct =>
        if ((ga == 0) || (gb == 0)) new BasisBlade()
        else if (Math.abs(ga - gb) == r.grade()) r
        else new BasisBlade()
      case ModifiedHestenesInnerProduct =>
        if (Math.abs(ga - gb) == r.grade()) r
        else new BasisBlade()
      case _ => null
    }

  def compareBlades(b1: BasisBlade, b2: BasisBlade): Int =
    if (b1.bitmap < b2.bitmap) -1
    else if (b1.bitmap > b2.bitmap) 1
    else java.lang.Double.compare(b1.scale, b2.scale)

  /** simplifies an arraylist (sum) of BasisBlades (destroys ordering of A!) */
  protected[clifford] def simplify(a: Seq[BasisBlade])(
    implicit m: Option[Metric] = None): Seq[BasisBlade] =
    if (a.size == 0) a
    else {
      val as = a.sortWith((x, y) => compareBlades(x, y) < 0)
      val result = Buffer[BasisBlade]()
      var current: BasisBlade = a(0).copy()
      for (b <- as) {
        if (b.bitmap == current.bitmap) {
          current = BasisBlade(current.bitmap, current.scale + b.scale)
        } else {
          if (current.scale != 0.0) result += current
          current = b.copy()
        }
      }
      if (current.scale != 0.0) result += current
      result.toSeq
    }

  protected[clifford] def roundDouble(
    what: Double, multipleOf: Double, epsilon: Double): Double = {
    val a: Double = what / multipleOf
    val b: Double = Math.round(a) * multipleOf
    if (Math.abs((what - b)) <= epsilon) b else what
  }
}

/**
  * @arg bitmap specifies what basis vectors are present in this basis blade
  * @arg scale scale of the basis blade is represented by this double
  */
case class BasisBlade(bitmap: Int, scale: Double) {
  /** constructs an instance of a unit BasisBlade */
  def this(b: Int) = this(b, 1.0)

  /** constructs an instance of a scalar BasisBlade */
  def this(s: Double) = this(0, s)

  /** constructs an instance of a zero BasisBlade */
  def this() = this(0, 0.0)

  def ⋅(other: BasisBlade)(implicit m: Option[Metric] = None): BasisBlade = 
    innerProduct(other)
  //def ·(other: BasisBlade)(implicit m: Option[Metric] = None): BasisBlade =
    //innerProduct(other)
  def ip(other: BasisBlade, tpe: InnerProductType = LeftContraction)(
    implicit m: Option[Metric] = None): BasisBlade =
    innerProduct(other)

  def innerProduct(other: BasisBlade, tpe: InnerProductType = LeftContraction)(
    implicit m: Option[Metric] = None): BasisBlade =
    BasisBlade.innerProduct(this, other, tpe)

  def ^(other: BasisBlade): BasisBlade = outerProduct(other)
  def op(other: BasisBlade): BasisBlade = outerProduct(other)

  def outerProduct(other: BasisBlade): BasisBlade =
    BasisBlade.outerProduct(this, other)

  def *(other: BasisBlade)(implicit m: Option[Metric] = None): Multivector =
    geometricProduct(other)
  def gp(other: BasisBlade)(implicit m: Option[Metric] = None): Multivector =
    geometricProduct(other)

  def geometricProduct(other: BasisBlade)(
    implicit m: Option[Metric] = None): Multivector = m match {
    case None => Multivector(BasisBlade.geometricProduct(this, other))
    case Some(metric) => {
      implicit val m1 = metric
      Multivector(BasisBlade.geometricProductM(this, other): _*)
    }
  }

  /**
    * @return reverse of this basis blade (always a newly constructed blade)
    */
  def reverse()(implicit m: Option[Metric] = None): BasisBlade =
    // multiplier = (-1)^(x(x-1)/2)
    BasisBlade(bitmap, minusOnePow((grade() * (grade() - 1)) / 2) * scale)

  // *!*HTML_TAG*!* grade_inversion
  /**
    * @return grade inversion of this basis blade (always a newly constructed 
    * blade)
    */
  def gradeInversion()(implicit m: Option[Metric] = None): BasisBlade =
    BasisBlade(bitmap, minusOnePow(grade()) * scale) // multiplier is (-1)^x

  // *!*HTML_TAG*!* clifford_conjugate
  /**
    * @return clifford conjugate of this basis blade (always a newly 
    * constructed blade)
    */
  def cliffordConjugate()(implicit m: Option[Metric] = None): BasisBlade =
    // multiplier is ((-1)^(x(x+1)/2)
    BasisBlade(bitmap, minusOnePow((grade() * (grade() + 1)) / 2) * scale)

  /** returns the grade of this blade */
  def grade(): Int = BasisBlade.bitCount(bitmap)

  // *!*HTML_TAG*!* minus_one_pow
  /** @return pow(-1, i) */
  def minusOnePow(i: Int): Int = if ((i & 1) == 0) 1 else -1

  //override def equals(o: Any): Boolean = o match {
    //case b: BasisBlade => ((b.bitmap == bitmap) && (b.scale == scale))
    //case _ => false
  //}

  //override def hashCode(): Int = java.util.Objects.hashCode(bitmap, scale)

  override def toString(): String = toString(null)

  /**
    * @param bvNames The names of the basis vector (e1, e2, e3) are used when
    * not available
    */
  def toString(bvNames: Array[String]): String = {
    val result: StringBuffer = new StringBuffer()
    var i = 1
    var b = bitmap
    while (b != 0) {
      if ((b & 1) != 0) {
	if (result.length() > 0) result.append("^")
	if (bvNames == null || i > bvNames.size || bvNames(i - 1) == null) {
	  result.append("e")
          result.append(i)
	} else {
          result.append(bvNames(i - 1))
        }
      }
      b >>= 1
      i = i + 1
    }

    if (result.length() == 0) scale.toString()
    else scale + "*" + result.toString()
  }

  /**
    * Rounds the scalar part of <code>this</code> to the nearest multiple X 
    * of <code>multipleOf</code>,
    * if |X - what| <= epsilon. This is useful when eigenbasis is used to 
    * perform products in arbitrary
    * metric, which leads to small roundof errors. You don't want to keep 
    * these roundof errors if your
    * are computing a multiplication table.
    *
    * @returns a new basis blade if a change is required.
    */
  def round(multipleOf: Double, epsilon: Double): BasisBlade = {
    val a: Double = BasisBlade.roundDouble(scale, multipleOf, epsilon)
    if (a != scale) BasisBlade(bitmap, a)
    else this
  }
}



