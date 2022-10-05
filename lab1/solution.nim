# sequence.sorted()
import std/algorithm

# ln(), sequence.sum(), sqrt(), exp()
import std/math

# initRand(), rand()
import std/random

# sequence.mapIt()
import std/sequtils

#############################
# Types
#############################

type
  StatisticValues = tuple[average, variance: float]
  HistogramUnit = tuple[lowerBound, upperBound: float, occurencies: int]
  Histogram = seq[HistogramUnit]

#############################
# Constants, global variables
#############################

const AMOUNT_OF_GENERATED_NUMBERS = 10000
const NUMBER_OF_FRACTIONAL_DIGITS = 2
const MINIMAL_NUMBER_OF_OCCURENCIES = 5

# For the current task, there is only one parameter of distribution law:
# lambda
const NUMBER_OF_DISTRIBUTION_LAW_PARAMETERS = 1

# Recommended amount of intervals
const AMOUNT_OF_INTERVALS = 20

# Level of significance equals 0.05
const PEARSONS_CHI_SQUARE_TEST = [
  3.8,
  6.0,
  7.8,
  9.5,
  11.1,
  12.6,
  14.1,
  15.5,
  16.9,
  18.3,
  19.7,
  21.0,
  22.4,
  23.7,
  25.0,
  26.3,
  27.6,
  28.9,
  30.1,
  31.4
]
var randomizer = initRand()

#############################
# Functions
#############################

proc randomNumber(): float = randomizer.rand(1.0)

proc givenFunction(lambd: float): float = (-1.0 / lambd) * ln(randomNumber())

proc calculateError(expected, actual: float): float =
  abs((expected - actual) / actual)

proc generateSequence(lambd: float): seq[float] =
  for _ in 1..AMOUNT_OF_GENERATED_NUMBERS:
    let number = givenFunction(lambd)
    result.add(number)


proc calculateStatisticValues(sequence: seq[float]): StatisticValues =
  let sumOfNumbers = sequence.sum()
  let amount = float(AMOUNT_OF_GENERATED_NUMBERS)

  let average = sumOfNumbers / amount

  let variancePart = sequence.mapIt((it - average)^2).sum()
  let squaredVariance = variancePart / (amount - 1.0)
  let variance = sqrt(squaredVariance)

  return (average, variance)


proc calculateHistogram(sequence: seq[float], intervalAmount: int): Histogram =
  let sortedSeq = sequence.sorted()

  let minimum = sortedSeq[sortedSeq.low]
  let maximum = sortedSeq[sortedSeq.high]
  let stepValue = (maximum - minimum) / float(intervalAmount)

  var intervalStep = minimum
  var index = sortedSeq.low

  for _ in 1..intervalAmount:
    let lowerBound = intervalStep
    intervalStep += stepValue
    let upperBound = intervalStep

    var occurencies = 0

    while true:
      var currentValue = sortedSeq[index]

      if currentValue <= intervalStep:
        inc occurencies
        inc index
        if index > sortedSeq.high: break
      else: break

    result.add((lowerBound: lowerBound, upperBound: upperBound, occurencies: occurencies))

proc joinIntervalsHelper (target: var HistogramUnit, source: HistogramUnit) =
  let intervalBounds = @[
    target.lowerBound,
    target.upperBound,
    source.lowerBound,
    source.upperBound
  ]

  let lowerBound = intervalBounds.min()
  let upperBound = intervalBounds.max()

  target.lowerBound = lowerBound
  target.upperBound = upperBound
  target.occurencies += source.occurencies

proc joinIntervals(histogram: var Histogram) =
  # Join intervals while any of them has a number of occurencies
  # lower than it's allowed
  while histogram.anyIt(it.occurencies < MINIMAL_NUMBER_OF_OCCURENCIES):
    var index = histogram.high

    while index > histogram.low:
      let currentInterval = histogram[index]
      let occurencies = currentInterval.occurencies

      if occurencies < MINIMAL_NUMBER_OF_OCCURENCIES:
        var intervalIndex: int

        # Join the interval with one which has lower number of occurencies
        if index == histogram.high:
          intervalIndex = index - 1
        else:
          let previousIndex = index - 1
          let nextIndex = index + 1
          intervalIndex = if histogram[previousIndex] > histogram[nextIndex]:
            nextIndex
            else: previousIndex

        joinIntervalsHelper(histogram[intervalIndex], currentInterval)
        histogram.delete(index)
        break
      dec index


proc hypoteticalDistributionLaw(x, lambd: float): float = 1 - exp(-lambd * x)

proc chiSquaredTest(histogram: Histogram, lambd: float): float =
  for interval in histogram:
    let lowerBound = interval.lowerBound
    let upperBound = interval.upperBound

    let theoreticalHitRate =
      hypoteticalDistributionLaw(upperBound, lambd) - hypoteticalDistributionLaw(lowerBound, lambd)
    let theoreticalOccurencies = (theoreticalHitRate * AMOUNT_OF_GENERATED_NUMBERS).round().int()

    result += (interval.occurencies - theoreticalOccurencies)^2 / theoreticalOccurencies



#############################
# Usage
#############################
const lambdas = [0.05, 0.1, 0.2, 0.25]

for lambd in lambdas:
  echo "\n============================="
  echo "Lambda: ", lambd

  let expectedAverage = 1.0 / lambd
  let expectedVariance = 1.0 / lambd

  let sequence = generateSequence(lambd)
  let statisticValues = calculateStatisticValues(sequence)

  echo "Expected average: ", expectedAverage
  echo "Actual average: ", statisticValues.average
  echo "Relative calculation error: ",
    calculateError(expectedAverage, statisticValues.average)

  echo "Expected variance: ", expectedVariance
  echo "Actual variance: ", statisticValues.variance
  echo "Relative calculation error: ",
    calculateError(expectedVariance, statisticValues.variance)

  var histogram = calculateHistogram(sequence, AMOUNT_OF_INTERVALS)
  histogram.joinIntervals()

  echo "Histogram sequence"
  for interval in histogram: echo interval.occurencies

  let chiSquaredTestValue = histogram.chiSquaredTest(lambd)
  echo "Chi-squared test: ", chiSquaredTestValue
  let freedomDegree = histogram.len - 1 - NUMBER_OF_DISTRIBUTION_LAW_PARAMETERS
  let theoreticalChiSquare = PEARSONS_CHI_SQUARE_TEST[freedomDegree]

  echo "Theoretical value of chi-square test: ", theoreticalChiSquare
  echo "Passed: ", chiSquaredTestValue < theoreticalChiSquare
