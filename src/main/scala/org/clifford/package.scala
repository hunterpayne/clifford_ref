
package org

package object clifford {

  implicit def arrayToMetric(arr: Array[Array[Double]]): Metric = Metric(arr)
}
