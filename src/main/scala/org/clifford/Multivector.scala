
package org.clifford

import scala.collection.mutable.Buffer

import breeze.linalg._

object Multivector {

  /** 
    * returns in the index [-1 31] of the highest bit that is on in 
    * <code>i</code> (-1 is returned if no bit is on at all (i == 0)) 
    */
  def highestOneBit(i: Int): Int = 31 - numberOfLeadingZeroBits(i)

  /** 
    * returns in the index [0 32] of the lowest bit that is on in 
    * <code>i</code> (32 is returned if no bit is on at all (i == 0))
    */
  def lowestOneBit(i: Int): Int = numberOfTrailingZeroBits(i)

  /** 
    * returns the number of 0 bits before the first 1 bit in <code>i</code>
    * <p>For example if i = 4 (100 binary), then 29 (31 - 2) is returned.
    */
  def numberOfLeadingZeroBits(arg: Int): Int = {
    // Note that unsigned shifting (>>>) is not required.
    var i = arg
    i = i | (i >> 1)
    i = i | (i >> 2)
    i = i | (i >> 4)
    i = i | (i >> 8)
    i = i | (i >>16)
    BasisBlade.bitCount(~i)
  }

  /**
    * returns the number of 0 bits after the last 1 bit in <code>i</code>
    * <p>For example if i = 4 (100 binary), then 2 is returned.
    */
  def numberOfTrailingZeroBits(i: Int): Int = BasisBlade.bitCount(~i & (i - 1))

  /** @return basis vector 'idx' range [0 ... dim) */
  def getBasisVector(idx: Int) = new Multivector(new BasisBlade(1 << idx))

  /** @return 'dim'-dimensional random vector with coordinates in 
    * range [-scale, scale] */
  def getRandomVector(dim: Int, scale: Double)(
    implicit m: Option[Metric] = None) =
    Multivector((for (i <- 0 until dim) yield {
      BasisBlade(1 << i, 2.0 * scale * (Math.random() - 0.5))
    }): _*)
    
  /**
    * @return 'dim'-dimensional random blade with coordinates in 
    * range [-scale, scale]
    */
  def getRandomBlade(dim: Int, grade: Int, scale: Double)(
    implicit m: Option[Metric] = None): Multivector = {
    var result = new Multivector(2 * scale * (Math.random() - 0.5))
    for (i <- 1 to grade) result = result.op(getRandomVector(dim, scale))
    result
  }

  def getRandomVersor(dim: Int, grade: Int, scale: Double)(
    implicit m: Option[Metric] = None) = {
    var result = new Multivector(2 * scale * (Math.random() - 0.5))
    for (i <- 1 to grade) result = result.gp(getRandomVector(dim, scale))
    result
  }

  protected def addToMatrix(
    m: DenseMatrix[Double], alpha: BasisBlade, beta: BasisBlade, 
    gamma: BasisBlade): Unit = 
    // add gamma.scale to matrix entry m[gamma.bitmap, beta.bitmap]
    m(gamma.bitmap, beta.bitmap) = m(gamma.bitmap, beta.bitmap) + gamma.scale

  protected def addToMatrix(
    m: DenseMatrix[Double], alpha: BasisBlade, beta: BasisBlade, 
    gamma: Seq[BasisBlade]): Unit =
    gamma.foreach { addToMatrix(m, alpha, beta, _) }

  def simplify(l: Seq[BasisBlade]): Seq[BasisBlade] = {
    val sorted = l.sortWith((x, y) => BasisBlade.compareBlades(x, y) < 0)
    if (sorted.isEmpty) Seq.empty
    else {
      val buf = Buffer[BasisBlade]()
      var prevBlade: BasisBlade = null
      sorted.foreach { b =>
        if (b.scale == 0.0) {
          if (prevBlade != null) buf += prevBlade
          prevBlade = null
        } else if (prevBlade != null && prevBlade.bitmap == b.bitmap) {
          prevBlade = new BasisBlade(b.bitmap, b.scale + prevBlade.scale)
        } else {
          if (prevBlade != null) buf += prevBlade
          prevBlade = b
        }
      }
      if (prevBlade != null) buf += prevBlade

      buf.filter { _.scale != 0.0 }.toSeq
    }
  }

  private def homogeneousPart(x: Multivector): (Multivector, Int) = {
    val grade = x.grade()
    if (grade < 0) {
      // homogeneousPart(x.largestGradePart())
      val newX = x.largestGradePart()
      (newX, newX.grade())
    } else {
      (x, grade)
    }
  }

  private def computeHomogeneousParts(
    aArg: Multivector, bArg: Multivector): 
      (Multivector, Int, Multivector, Int) = {

    val (a, ga) = homogeneousPart(aArg)
    val (b, gb) = homogeneousPart(bArg)

    // normalize (approximately) and swap (optionally)
    if (ga <= gb) {
      (a.gpD(1.0 / a.largestCoordinate().toDouble), ga, 
        b.gpD(1.0 / b.largestCoordinate().toDouble), gb)
    } else {
      (b.gpD(1.0 / b.largestCoordinate().toDouble), gb,
        a.gpD(1.0 / a.largestCoordinate().toDouble), ga)
    }
  }

  def meet(a: Multivector, b: Multivector): Multivector = {

    val largeEpsilon: Double = 10e-4
    val (ca, ga, cb, gb) = computeHomogeneousParts(a, b)

    // compute delta product & 'normalize'
    val dScore = deltaProduct(ca, cb)
    val gd = dScore.grade()
    val d = dScore.gpD(1.0 / dScore.largestCoordinate().toDouble)

    // if delta product is scalar, we're done:
    if (0 == gd) {
      ca
      // if grade of delta product is equal to ga + gb, we're done, too
    } else if (gd == ga + gb) {
      new Multivector(1.0)
    } else {
      val dim = Math.max(ca.spaceDim(), cb.spaceDim())
      // the pseudoscalar
      val i = new Multivector(new BasisBlade((1 << dim) - 1, 1.0))

      // init join to pseudoscalar
      var j = i
      // compute excessity of join
      var ej = dim - ((ga + gb + gd) >> 1)

      var m = new Multivector(1.0)
      var em = ((ga + gb - gd) >> 1) // compute excessity of meet

      // init s, the dual of the delta product:
      val s: Multivector = d.ip(i.versorInverse(), LeftContraction)

      // precompute inverse of ca
      val cai: Multivector = ca.versorInverse()
      for (i <- 0 until dim) {
	val c = (new Multivector(new BasisBlade(1 << i, 1.0))).
          ip(s, ModifiedHestenesInnerProduct).
          ip(s, ModifiedHestenesInnerProduct)

        if (c.largestCoordinate() >= largeEpsilon) {
	  // compute projection, rejection of 'c' wrt to 'ca'
          // use correct inverse because otherwise cr != c - cp
          val cp = c.ip(ca, LeftContraction).ip(cai, LeftContraction)
          val cr = c - cp

	  // if 'c' has enough of it in 'ca', then add to meet
	  if (cp.largestCoordinate() > largeEpsilon) {
	    m = m.op(cp)
	    em = em - 1 // reduce excessity of meet
	    if (em == 0) { // has the meet been computed???
	      return m
	    }

	    if (cr.largestCoordinate() > largeEpsilon) {
	      j = cr.ip(j, LeftContraction)
	      ej = ej - 1 // reduce excessity of join
	      if (0 == ej) { // has the join been computed???
                return cb.ip(j.versorInverse(), LeftContraction).
                  ip(ca, LeftContraction)
              }
            }
	  }
        }
      }

      throw new java.lang.ArithmeticException("meet algorithm failed!")
    }
  }

  def join(a: Multivector, b: Multivector): Multivector = {
    val largeEpsilon: Double = 10e-4
    val (ca, ga, cb, gb) = computeHomogeneousParts(a, b)

    // compute delta product & 'normalize'
    val dScore = deltaProduct(ca, cb)
    val gd = dScore.grade()
    val d = dScore.gpD(1.0 / dScore.largestCoordinate().toDouble)

    // if delta product is scalar, we're done:
    if (0 == gd) {
      ca ^ ca.versorInverse().ip(cb, LeftContraction)
      // if grade of delta product is equal to ga + gb, we're done, too
    } else if (gd == ga + gb) {
      ca ^ cb
    } else {
      val dim = Math.max(ca.spaceDim(), cb.spaceDim())
      // the pseudoscalar
      val i = new Multivector(new BasisBlade((1 << dim) - 1, 1.0))

      // init join to pseudoscalar
      var j = i
      // compute excessity of join
      var ej = dim - ((ga + gb + gd) >> 1)

      var m = new Multivector(1.0)
      var em = ((ga + gb - gd) >> 1) // compute excessity of meet

      // init s, the dual of the delta product:
      val s: Multivector = d.ip(i.versorInverse(), LeftContraction)

      // precompute inverse of ca
      val cai: Multivector = ca.versorInverse()
      for (i <- 0 until dim) {
	val c = (new Multivector(new BasisBlade(1 << i, 1.0))).
          ip(s, ModifiedHestenesInnerProduct).
          ip(s, ModifiedHestenesInnerProduct)

        if (c.largestCoordinate() >= largeEpsilon) {
	  // compute projection, rejection of 'c' wrt to 'ca'
          // use correct inverse because otherwise cr != c - cp
          val cp = c.ip(ca, LeftContraction).ip(cai, LeftContraction)
          val cr = c - cp

	  // if 'c' has enough of it in 'ca', then add to meet
	  if (cp.largestCoordinate() > largeEpsilon) {
	    m = m ^ cp
	    em = em - 1 // reduce excessity of meet
	    if (em == 0) { // has the meet been computed???
	      return ca ^ ca.versorInverse().ip(cb, LeftContraction)
	    }

	    if (cr.largestCoordinate() > largeEpsilon) {
	      j = cr.ip(j, LeftContraction)
	      ej = ej - 1 // reduce excessity of join
	      if (0 == ej) { // has the join been computed???
                return j
              }
            }
	  }
        }
      }

      throw new java.lang.ArithmeticException("meet algorithm failed!")
    }
  }

  def deltaProduct(a: Multivector, b: Multivector)(
    implicit m: Option[Metric] = None): Multivector = {
    val d = (a * b).compress()
    d.extractGrade(d.topGradeIndex())
  }

  def apply(blades: BasisBlade*): Multivector = 
    new Multivector(
      simplify(blades).
        sortWith((x, y) => BasisBlade.compareBlades(x, y) < 0): _*)
}

class Multivector(val blades: BasisBlade*) {

  import Multivector._

  def this(s: Double) = this(new BasisBlade(s))

  override def equals(o: Any): Boolean = o match {
    case m: Multivector => {
      val zero = subtract(m)
      zero.blades.size == 0
    }
    case _ => false
  }

  override def toString: String = toString(null)

  /**
    * @param bvNames The names of the basis vector (e1, e2, e3) are used when
    * not available
    */
  def toString(bvNames: Array[String]): String =
    if (blades.size == 0) {
      "0"
    } else {
      val result = new StringBuffer()
      for ((b, i) <- blades.zipWithIndex) {
	val s = b.toString(bvNames)
	if (i == 0) {
          result.append(s)
	} else if (s(0) == '-') {
	  result.append(" - ")
	  result.append(s.substring(1))
	} else {
	  result.append(" + ")
	  result.append(s)
	}
      }
      result.toString()
    }

  /** @return geometric product of this with a scalar */
  def gpD(a: Double)(implicit m: Option[Metric] = None): Multivector =
    if (a == 0.0) {
      new Multivector()
    } else {
      Multivector(
        (for (b <- blades) yield new BasisBlade(b.bitmap, b.scale * a)): _*)
    }

  def *(x: Multivector)(implicit m: Option[Metric] = None): Multivector = gp(x)

  // *!*HTML_TAG*!* gp
  /** @return geometric product of this with a 'm' using metric 'M' */
  def gp(x: Multivector)(implicit m: Option[Metric] = None): Multivector = {
    val newBlades = for (b1 <- blades; b2 <- x.blades) yield((b1 * b2).blades)
    Multivector(newBlades.flatten: _*)
  }

  def ^(x: Multivector)(implicit m: Option[Metric] = None): Multivector = op(x)
   
  // *!*HTML_TAG*!* op
  /** @return outer product of this with 'x' */
  def op(x: Multivector)(implicit m: Option[Metric] = None): Multivector = 
    Multivector((for (b1 <- blades; b2 <- x.blades) yield (b1 ^ b2)): _*)

  def ⋅(other: Multivector)(implicit m: Option[Metric] = None): Multivector =
    ip(other)

  def x(other: Multivector)(implicit m: Option[Metric] = None): Multivector =
    ip(other)

  def ip(x: Multivector, tpe: InnerProductType = LeftContraction)(
    implicit m: Option[Metric] = None): Multivector = 
    Multivector((for (b1 <- blades; b2 <- x.blades) yield b1.ip(b2, tpe)): _*)

  /** @return scalar product of this with a 'x' */
  def scalarProduct(x: Multivector)(implicit m: Option[Metric] = None): Double =
    ip(x, LeftContraction).scalarPart()

  /** shortcut to scalarProduct(...) */
  def scp(x: Multivector)(implicit m: Option[Metric] = None): Double = 
    scalarProduct(x)

  //def +(a: Double)(implicit m: Option[Metric] = None): Multivector = addD(a)
  def +(x: Multivector)(implicit m: Option[Metric] = None): Multivector = add(x)
  //def -(a: Double)(implicit m: Option[Metric] = None): Multivector = 
    //subtractD(a)
  def -(x: Multivector)(implicit m: Option[Metric] = None): Multivector = 
    subtract(x)

  /** @return sum of this with scalar 'a' */
  def addD(a: Double)(implicit m: Option[Metric] = None): Multivector =
    Multivector((blades.toSeq :+ (new BasisBlade(a))).map { _.copy() }: _*)

  /** @return sum of this with 'x' */
  def add(x: Multivector)(implicit m: Option[Metric] = None): Multivector =
    Multivector((blades ++ x.blades).map { _.copy() }: _*)

  /** @return this - scalar 'a' */
  def subtractD(a: Double)(implicit m: Option[Metric] = None): Multivector =
    addD(-a)

  /** @return this - 'x' */
  def subtract(x: Multivector)(implicit m: Option[Metric] = None): Multivector =
    Multivector((blades ++ x.gpD(-1).blades).map { _.copy() }: _*)

  def ∪(x: Multivector)(implicit m: Option[Metric] = None): Multivector = 
    join(x)

  def join(x: Multivector)(implicit m: Option[Metric] = None): Multivector =
    Multivector.join(this, x)

  def ∩(x: Multivector)(implicit m: Option[Metric] = None): Multivector = 
    meet(x)

  def meet(x: Multivector)(implicit m: Option[Metric] = None): Multivector =
    Multivector.meet(this, x)

  // *!*HTML_TAG*!* exp
  /** @return exponential of this */
  def exp(order: Int = 12)(implicit m: Option[Metric] = None): Multivector = {
    // check out this^2 for special cases
    val a2: Multivector = (this * this).compress()
    if (a2.isNull(1e-8)) {
      // special case A^2 = 0
      addD(1)
    } else if (a2.isScalar()) {
      val a2s: Double = a2.scalarPart()
      // special case A^2 = +-alpha^2
      if (a2s < 0) {
	val alpha = Math.sqrt(-a2s)
	this.gpD(Math.sin(alpha) / alpha).addD(Math.cos(alpha))
      } else {
	val alpha: Double = Math.sqrt(a2s)
	this.gpD(MathU.sinh(alpha) / alpha).addD(MathU.cosh(alpha))
      } // TODO what if a2 == 0?
    } else {
      expSeries(order)
    }
  }

  // *!*HTML_TAG*!* expSeries
  /** Evaluates exp using series . . .  (== SLOW & INPRECISE!) */
  protected def expSeries(order: Int)(
    implicit m: Option[Metric] = None): Multivector = {
    // first scale by power of 2 so that its norm is ~ 1
    var scale: Long = 1
    var max: Double = norm_e()
    while (max > 1.0) {
      max = max / 2.0
      scale = scale << 1
    }

    val scaled: Multivector = gpD(1.0 / scale.toDouble)

    // taylor approximation
    var result = new Multivector(1.0)
    var tmp = new Multivector(1.0)

    for (i <- 1 until order) {
      tmp = tmp.gp(scaled.gpD(1.0 / i.toDouble))
      result = Multivector((result.blades ++ tmp.blades): _*)
    }

    // undo scaling
    while (scale > 1) {
      result = result * result
      scale = scale >>> 1
    }
    result
  }

  // *!*HTML_TAG*!* sin
  def sin(order: Int = 12)(implicit m: Option[Metric] = None): Multivector = {
    // check out this^2 for special cases
    val a2: Multivector = (this * this).compress()
    if (a2.isNull(1e-8)) {
      // special case A^2 = 0
      this
    } else if (a2.isScalar()) {
      val a2s: Double = a2.scalarPart()
      // special case A^2 = +-alpha^2
      if (a2s < 0.0) {
        val alpha = Math.sqrt(-a2s)
        gpD(MathU.sinh(alpha) / alpha)
      } else {
        val alpha = Math.sqrt(a2s)
        gpD(Math.sin(alpha) / alpha)
      }
    } else {
      sinSeries(order)
    }
  }

  // *!*HTML_TAG*!* sinSeries
  /** Evaluates sin using series . . .  (== SLOW & INPRECISE!) */
  protected def sinSeries(order: Int)(
    implicit m: Option[Metric] = None): Multivector = {

    var result = this
    var tmp = this
    var sign: Int = -1

    for (i <- 2 until order) {
      tmp = tmp * this.gpD(1.0 / i.toDouble)
      if ((i & 1) != 0) { // only odds
        result = Multivector(
          (result.blades ++ (tmp.gpD(sign.toDouble)).blades): _*)
        sign = -sign
      }
    }

    result
  }

  // *!*HTML_TAG*!* cos
  def cos(order: Int = 12)(implicit m: Option[Metric] = None): Multivector = {
    // check out this^2 for special cases
    val a2: Multivector = (this * this).compress()
    if (a2.isNull(1e-8)) {
      // special case A^2 = 0
      new Multivector(1)
    } else if (a2.isScalar()) {
      val a2s: Double = a2.scalarPart()
      // special case A^2 = +-alpha^2
      if (a2s < 0.0) {
        val alpha = Math.sqrt(-a2s)
        new Multivector(MathU.cosh(alpha))
      } else {
        val alpha = Math.sqrt(a2s)
        new Multivector(Math.cos(alpha))
      }
    } else {
      cosSeries(order)
    }
  }

  // *!*HTML_TAG*!* cosSeries
  /** Evaluates cos using series . . .  (== SLOW & INPRECISE!) */
  protected def cosSeries(order: Int)(
    implicit m: Option[Metric] = None): Multivector = {

    var result = new Multivector(1.0)
    var tmp = this
    var sign: Int = -1

    for (i <- 2 until order) {
      tmp = tmp * this.gpD(1.0 / i.toDouble)
      if ((i & 1) == 0) { // only the even part of the series
        result = Multivector(
          (result.blades ++ tmp.gpD(sign.toDouble).blades): _*)
        sign = -sign
      }
    }

    result
  }


  /**
    * Can throw java.lang.ArithmeticException if multivector is null
    * @return unit under Euclidean norm
    */
  def unit_e(): Multivector = unit_r()

  def norm_e(): Double = {
    val s = scp(reverse())
    if (s < 0.0) 0.0
    else Math.sqrt(s)
  }

  def norm_e2(): Double = math.max(0.0, scp(reverse()))

  /**
    * Can throw java.lang.ArithmeticException if multivector is null
    * @return unit under 'reverse' norm (this / sqrt(abs(this.reverse(this))))
    */
  def unit_r()(implicit m: Option[Metric] = None): Multivector = {
    val s = scp(reverse())
    if (0.0 == s) throw new Exception("null multivector")
    else gpD(1.0 / Math.sqrt(Math.abs(s)))
  }

  //public Multivector unit_r() {
    //double s = scp(reverse());
    //if (s == 0.0) throw new java.lang.ArithmeticException("null multivector");
    //else return this.gp(1 / Math.sqrt(Math.abs(s)));
  //}

  /** @return true if this is really 0.0 */
  def isNull(): Boolean = {
    simplifyMe()
    blades.isEmpty
  }

  /** @return true if norm_e2 < epsilon * epsilon*/
  def isNull(epsilon: Double): Boolean = {
    val s = norm_e2()
    (s < epsilon * epsilon)
  }

  /** @return true is this is a scalar (0.0 is also a scalar) */
  def isScalar(): Boolean = 
    if (isNull()) true
    else if (blades.size == 1) blades(0).bitmap == 0
    else false

  // *!*HTML_TAG*!* reverse
  /** @return reverse of this */
  def reverse(): Multivector = Multivector((blades.map { _.reverse() }): _*)

  // *!*HTML_TAG*!* grade_inversion
  /** @return grade inversion of this */
  def gradeInversion(): Multivector = 
    Multivector((blades.map { _.gradeInversion() }): _*)

  // *!*HTML_TAG*!* clifford_conjugate
  /** @return clifford conjugate of this */
  def cliffordConjugate(): Multivector =
    Multivector((blades.map { _.cliffordConjugate() }): _*)

  /**
    * Extracts grade 'g' from this multivector.
    * @return a new multivector of grade 'g'
    */
  def extractGrade(grades: Int*): Multivector = {
    // what is the maximum grade to be extracted?
    require(grades.size > 0)
    Multivector((blades.filter { b => grades.contains(b.grade()) }): _*)
  }

  def dual(dim: Int)(implicit metric: Option[Metric] = None): Multivector = {
    val i = new Multivector(new BasisBlade((1 << dim) - 1, 1.0))
    ip(i.versorInverse(), LeftContraction)
  }

  def dual()(implicit metric: Option[Metric]): Multivector = {
    // can't use dual void arg variant without a metric
    require(metric.isDefined)
    val i = new Multivector(
      new BasisBlade((1 << metric.get.eigenMetric.size) - 1, 1.0))
    ip(i.versorInverse(), LeftContraction)
  }

  def scalarPart()(implicit m: Option[Metric] = None): Double = 
    blades.filter { _.bitmap == 0 }.map { _.scale }.sum

  // *!*HTML_TAG*!* versor_inverse
  /**
    * Can throw java.lang.ArithmeticException if versor is not invertible
    * @return inverse of this (assuming it is a versor, no check is made!)
    */
  def versorInverse()(implicit m: Option[Metric] = None): Multivector = {
    val r: Multivector = reverse()
    val s: Double = scp(r)
    if (s == 0.0) 
      throw new java.lang.ArithmeticException("non-invertible multivector")
    r.gpD(1.0 / s)
  }

  // *!*HTML_TAG*!* general_inverse
  /**
    * Can throw java.lang.ArithmeticException if blade is not invertible
    * @return inverse of arbitrary multivector.
    *
    */
  def generalInverse()(implicit m: Option[Metric] = None): Multivector = {
    val dim = spaceDim()
    val width = 1 << dim

    val matrix: DenseMatrix[Double] = DenseMatrix.zeros(width, width)

    // create all unit basis blades for 'dim'
    val b = (0 until width).map { new BasisBlade(_) }

    // construct a matrix 'M' such that matrix multiplication of 'M' with
    // the coordinates of another multivector 'x' (stored in a vector)
    // would result in the geometric product of 'M' and 'x'
    blades.foreach { blade => (0 until width).foreach { j =>
      addToMatrix(matrix, blade, b(j), (blade * b(j)).blades)
    }}

    // try to invert matrix (can fail, then we throw an exception)
    val im = inv(matrix)

    // reconstruct multivector from first column of matrix
    Multivector(
      ((0 until width).map { j => (j, im(j, 0)) }.
        filter { case(j, v) => v != 0.0 }.
        map { case(j, v) => b(j).copy(scale = v) }): _*)
  }

  /** @return simplification of this multivector (the same Multivector,
    * but blades array can be changed) */
  def simplifyMe(): Multivector = {
    //simplify(blades)
    this
  }

  /** @return abs(largest coordinate) of 'this' */
  def largestCoordinate(): Double = {
    //simplifyMe()
    if (!blades.isEmpty) blades.map { b => Math.abs(b.scale) }.max
    else 0
  }

  /** @return abs(largest BasisBlade) of 'this' */
  def largestBasisBlade(): BasisBlade =
    blades.map { b => (b, Math.abs(b.scale)) }.maxBy { _._2 }._1

  /** 
    * @return the grade of this if homogeneous, -1 otherwise.
    * 0 is return for null Multivectors.
    */
  def grade(): Int = {
    val uniqGrades = blades.map { _.grade() }.distinct
    if (uniqGrades.isEmpty) 0
    else if (uniqGrades.size > 1) -1
    else uniqGrades(0)
  }

  /** @return bitmap of grades that are in use in 'this'*/
  def gradeUsage(): Int =
    blades.map { _.grade() }.foldLeft(0)(_ | 1 << _)

  /** @return index of highest grade in use in 'this'*/
  def topGradeIndex(): Int = blades.map { _.grade() }.max

  /** @return the largest grade part of this */
  def largestGradePart(): Multivector = {
    simplifyMe()

    var maxGP: Multivector = null
    var maxNorm: Double = -1.0
    val gu = gradeUsage()
    for { i <- 0 to topGradeIndex() } {
      if ((gu & (1 << i)) != 0) {
        val gp: Multivector = extractGrade(i)
        val n: Double = gp.norm_e()
        if (n > maxNorm) {
          maxGP = gp
          maxNorm = n
        }
      }
    }

    if (maxGP == null) new Multivector() else maxGP
  }

  /** @return dimension of space this blade (apparently) lives in */
  protected def spaceDim(): Int = 
    blades.map { b => highestOneBit(b.bitmap) }.max + 1

  /**
    * Currently removes all basis blades with |scale| less than epsilon
    *
    * Old version did this:
    * Removes basis blades with whose |scale| is less than 
    * <code>epsilon * maxMag</code> where
    * maxMag is the |scale| of the largest basis blade.
    *
    * @return 'Compressed' version of this (the same Multivector, but 
    * blades array can be changed)
    */
  def compress(epsilon: Double = 1e-13): Multivector = 
    if (largestCoordinate() == 0.0) {
      new Multivector()
    } else {
      // premultiply maxMag
      //maxMag = epsilon // used to read *=
      Multivector(blades.filter { b => Math.abs(b.scale) >= epsilon }: _*)
    }



  /** sorts the blade in 'blades' based on bitmap only */
  //protected def sortBlades(): Unit = {
    //if (bladesSorted) return;
    //else {
      //Collections.sort(blades, new BladesComperator());
      //bladesSorted = true;
    //}
  //}

  
}
