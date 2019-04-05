
package org.clifford

import org.scalatest.{ FlatSpec, Matchers }

class MultivectorSpace extends FlatSpec with Matchers {

  behavior of "Multivector"

  it should "bit operations should work" in {

    import Multivector._

    numberOfLeadingZeroBits(0x40000000) should be (1)
    numberOfLeadingZeroBits(4) should be (29)
    numberOfTrailingZeroBits(4) should be (2)
    numberOfTrailingZeroBits(0x40000000) should be (30)
    lowestOneBit(4) should be (2)
    lowestOneBit(6) should be (1)
    lowestOneBit(5) should be (0)
    lowestOneBit(8) should be (3)
    lowestOneBit(0) should be (32)
    highestOneBit(7) should be (2)
    highestOneBit(0x40000000) should be (30)
  }

  it should "complete basic operations" in {
  
    import Multivector._

    // setup conformal algebra:
    val bvNames: Array[String] = Array("no", "e1", "e2", "e3", "ni")
    val m: Array[Array[Double]] = Array(
      Array(0.0, 0.0, 0.0, 0.0, -1.0),
      Array(0.0, 1.0, 0.0, 0.0, 0.0),
      Array(0.0, 0.0, 1.0, 0.0, 0.0),
      Array(0.0, 0.0, 0.0, 1.0, 0.0),
      Array(-1.0, 0.0, 0.0, 0.0, 0.0)
    )

    implicit val metric: Metric = Metric(m)

    val no: Multivector = getBasisVector(0)
    val ni: Multivector = getBasisVector(4)
    val e1: Multivector = getBasisVector(1)
    val e2: Multivector = getBasisVector(2)
    val e3: Multivector = getBasisVector(3)

    val a: Multivector = e1.add(e2.op(e3).op(e1))
    println("A = " + a)
    //println(new MultivectorType(A));

    val a2 = e1.add(e2.op(e3))
    val ai1 = a2.generalInverse()
    val t1 = ai1.gp(a2).subtractD(1.0).compress(1e-7)
    val t2 = a2.gp(ai1).subtractD(1.0).compress(1e-7)
    //t1.isNull() should be (false)
    //t2.isNull() should be (false)

    val dim = 5
    for (i <- 0 until 1000) {

      val a = Multivector.getRandomVersor(dim, (Math.random() * dim).toInt, 1.0)
      val ai1 = a.generalInverse()
      val ai2 = a.versorInverse()

      //(ai1 ^ ai2).isNull should be (false)

      val t1 = ai1.gp(a).subtractD(1.0).compress(1e-7)
      val t2 = a.gp(ai1).subtractD(1.0).compress(1e-7)

      //t1.isNull() should be (false)
      //t2.isNull() should be (false)
    }

    val b = new Multivector(-1.45)
    val r = b.cos()
    //val r2 = b.cosSeries(24)

    println("B = " + b.toString(bvNames))
    println("R1 = " + r.toString(bvNames))
    //println("R2 = " + r2.toString(bvNames))

    val b2 = e1.op(e3).gpD(1.33334)
    val r3 = b2.cos()
    //val r4 = b2.cosSeries(24)

    println("B = " + b2.toString(bvNames))
    println("R1 = " + r3.toString(bvNames))
    //println("R2 = " + r4.toString(bvNames))
  }

  it should "complete meet operations" in {

    // test code for meet / join
    val dim = 7

    for (i <- 0 until 1000) {
      var a = Multivector.getRandomBlade(
        dim, (0.6 * Math.random() * dim).toInt, 1.0)
      var b = Multivector.getRandomBlade(
        dim, (0.6 * Math.random() * dim).toInt, 1.0)

      if (Math.random() < 0.25) {
	// make a subspace of b:
	b = a.op(Multivector.getRandomBlade(
          dim, (0.5 * Math.random() * dim + 0.99).toInt, 1.0))
      }

      if (Math.random() < 0.25) {
	// use basis elements for 'a'
	a = new Multivector(
          new BasisBlade((Math.random() * ((1 << 7) - 1)).toInt, 1.0))
      }
      if (Math.random() < 0.25) {
	// use basis elements for 'b'
	b = new Multivector(
          new BasisBlade((Math.random() * ((1 << 7) - 1)).toInt, 1.0))
      }

      if (Math.random() < 0.5) {
	// swap a & b:
	val tmp = b
	b = a
	a = tmp
      }

      if (!a.isNull(1e-3) && !b.isNull(1e-3)) {

        //	    println("a = " + a + ",");
        //	    println("b = " + b + ",");
        val met = a.meet(b)

        met.ip(a, LeftContraction).compress(1e-7).isNull() should be (false)
        met.ip(b, LeftContraction).compress(1e-7).isNull() should be (false)
      }

      println("Loop " + i + " ")
    }
  }

  it should "complete join operations" in {

    // test code for meet / join
    val dim = 7

    // TODO join is broken
    for (i <- 0 until 0) {
      var a = Multivector.getRandomBlade(
        dim, (0.6 * Math.random() * dim).toInt, 1.0)
      var b = Multivector.getRandomBlade(
        dim, (0.6 * Math.random() * dim).toInt, 1.0)

      if (Math.random() < 0.25) {
	// make a subspace of b:
	b = a.op(Multivector.getRandomBlade(
          dim, (0.5 * Math.random() * dim + 0.99).toInt, 1.0))
      }

      if (Math.random() < 0.25) {
	// use basis elements for 'a'
	a = new Multivector(
          new BasisBlade((Math.random() * ((1 << 7) - 1)).toInt, 1.0))
      }
      if (Math.random() < 0.25) {
	// use basis elements for 'b'
	b = new Multivector(
          new BasisBlade((Math.random() * ((1 << 7) - 1)).toInt, 1.0))
      }

      if (Math.random() < 0.5) {
	// swap a & b:
	val tmp = b
	b = a
	a = tmp
      }

      if (!a.isNull(1e-3) && !b.isNull(1e-3)) {

        //	    println("a = " + a + ",");
        //	    println("b = " + b + ",");
        val met = a.join(b)

        met.ip(a, LeftContraction).compress(1e-7).isNull() should be (false)
        met.ip(b, LeftContraction).compress(1e-7).isNull() should be (false)
      }

      println("Loop " + i + " ")
    }
  }

  it should "handle complex operations" in {

    /*
    // setup conformal algebra:
    //String[] bvNames = {"no", "e1", "e2", "e3", "e4", "ni"};
    double[][] m = new double[][]{
      {0.0, 0.0, 0.0, 0.0, -1.0},
      {0.0, 1.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 1.0, 0.0},
      {-1.0, 0.0, 0.0 , 0.0, 0.0}
    };

    Metric M = null;
    try {
      M = new Metric(m);
    } catch (Exception ex) {}


    Multivector no = Multivector.getBasisVector(0);
    Multivector e1 = Multivector.getBasisVector(1);
    Multivector e2 = Multivector.getBasisVector(2);
    Multivector e3 = Multivector.getBasisVector(3);
    //	Multivector e4 = Multivector.getBasisVector(4);
    Multivector ni = Multivector.getBasisVector(4);

    long t = System.currentTimeMillis();
    // test code for factorization
    int dim = 8;
    double[] scale = new double[1];
    ArrayList[] SSS = new ArrayList[dim+1];
    for (int i = 0; i <= dim; i++)
      SSS[i] = new ArrayList();

    for (int i = 0; i < 1; i++) {
      Multivector B = Multivector.getRandomBlade(dim, (int)(Math.random() * (dim + 0.49)), 1.0);

      ArrayList BL = new ArrayList();
      BL.add(new BasisBlade(30, -0.662244));
      BL.add(new BasisBlade(29, -0.391495));
      BL.add(new BasisBlade(27, -0.430912));
      BL.add(new BasisBlade(23, 0.218277));
      BL.add(new BasisBlade(15, -0.213881));
      B = new Multivector(BL);

      Multivector[] f = factorizeBlade(B, scale);
      Multivector R = new Multivector(scale[0]);
      for (int g = 0; g < f.length; g++)
	R = R.op(f[g]);

      Multivector[] fAltFast = factorizeBladeAltFast(B, scale);
      Multivector RaltFast = new Multivector(scale[0]);
      for (int g = 0; g < fAltFast.length; g++) {
        //                            System.out.println("f: " + fAltFast[g]);
        RaltFast = RaltFast.op(fAltFast[g]);
      }

      B = B.unit_e();
      R = R.unit_e();
      RaltFast = RaltFast.unit_e();

      int k = B.grade();

      double checkScale = R.gp(RaltFast.versorInverse()).scalarPart();
      if (checkScale < 0)
        System.out.println("Whaaaaa!\n");


      SSS[B.grade()].add(
        (checkScale < 0) ? "-" : "+");


      System.out.println("B  = " + B + ",");
      System.out.println("R  = " + R + ",");
      System.out.println("Ra = " + RaltFast + ",");
      /*
       if ((i % 100) == 0)
       System.out.println("I: " + i);*/

    }
    for (int i = 0; i <= dim; i++)
      System.out.println("Dim " + i + " = " + SSS[i].toString());
    System.out.println("Done!" + (System.currentTimeMillis() - t));
     */


    //Multivector V = e1.add(e2).gp(no.add(ni)).gp(e3);
    //V = e1.add(e2).gp(no.add(ni));
    /*for (int i = 0; i < 100; i++) {
     int dim = M.getEigenMetric().length;
     Multivector V = Multivector.getRandomVersor(dim, (int)(Math.random() * (dim + 0.49)), 1.0, M);
     //System.out.println("V = " + V.toString(bvNames));
     factorizeVersor(V, M);
     }*/


    /*
     double alpha = 0.1 * Math.PI / 2;
     double[][] m = new double[][]{
     new double[]{Math.cos(alpha), Math.sin(alpha), 0.0},
     new double[]{-Math.sin(alpha), Math.cos(alpha), 0.0},
     new double[]{0.0, 0.0, 1.0}
     };*/
    /*
     double[][] m = new double[][]{
     new double[]{1.0, 0.0, 0.0},
     new double[]{0.0, 1.0, 0.0},
     new double[]{0.0, 0.0, 1.0}
     };
     */
    /*
     Multivector R = rotationMatrixToRotor(m, e1, e2, e3);
     Multivector Ralt = rotationMatrixToRotorAlt(m, e1, e2, e3);
     
     System.out.println("R = " + R + ",");
     System.out.println("Ralt = " + Ralt + ",");
     */
  }
}

