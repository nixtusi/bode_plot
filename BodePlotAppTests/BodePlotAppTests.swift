//
//  BodePlotAppTests.swift
//  BodePlotAppTests
//
//  Created by Yuta Nisimatsu on 2025/12/23.
//

import Testing
import SwiftUI
import UIKit
@testable import BodePlotApp

struct BodePlotAppTests {

    @Test(arguments: ["1.5", "1,5", "١٫٥", "１．５"])
    func localizedDecimalInput(input: String) throws {
        var lexer = Lexer(input)
        #expect(try lexer.nextToken() == .number(1.5))
        #expect(try lexer.nextToken() == .end)
    }

    @Test func localizedDecimalInTransferFunction() throws {
        let expression = try Parser("١٫٥/(s+１)").parse()
        let coefficients = expression.toHighFirstCoeffs()
        #expect(coefficients.num == [1.5])
        #expect(coefficients.den == [1, 1])
    }

    @Test func calculationReflectsEditsAndSettings() async throws {
        let calculator = BodeCalculator()
        let first = ExprEntry(text: "1/(s+1)", colorIndex: 0)
        var second = ExprEntry(text: "sqrt(s+1)", colorIndex: 1)

        let initial = try #require(await calculator.calculate(entries: [first, second], start: -1, end: 1, points: 3))
        #expect(initial.count == 2)
        let pole = try #require(initial[0].bodeData.dropFirst().first)
        #expect(abs(pole.w - 1) < 1e-12)
        #expect(abs(pole.magDB - (-3.0102999566)) < 1e-8)
        #expect(abs(pole.phaseDeg - (-45)) < 1e-10)
        let squareRoot = try #require(initial[1].bodeData.dropFirst().first)
        #expect(abs(squareRoot.magDB - 1.5051499783) < 1e-8)
        #expect(abs(squareRoot.phaseDeg - 22.5) < 1e-10)

        second.text = "sqrt("
        let edited = try #require(await calculator.calculate(entries: [first, second], start: -1, end: 1, points: 3))
        #expect(edited[0].bodeData.map(\.magDB) == initial[0].bodeData.map(\.magDB))
        #expect(edited[1].bodeData.isEmpty)
        #expect(edited[1].errorMessage != nil)

        second.text = "sqrt(s+1)"
        let resized = try #require(await calculator.calculate(entries: [second], start: -2, end: 2, points: 5))
        #expect(resized[0].bodeData.count == 5)
        #expect(resized[0].bodeData.first?.w == 0.01)
        #expect(resized[0].bodeData.last?.w == 100)
        #expect(resized[0].label == "G₁")
        #expect(resized[0].bodeData.allSatisfy { $0.seriesLabel == "G₁" })
    }

    @Test func cancelledCalculationDoesNotPublishResults() async {
        let calculator = BodeCalculator()
        let task = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return await calculator.calculate(entries: [ExprEntry(text: "1/(s+1)")],
                                              start: -2, end: 2, points: 3000)
        }
        if await task.value != nil {
            Issue.record("A cancelled calculation must not produce results for the UI")
        }
    }

    @Test func largePowersAreEvaluatedWithoutExpandingHugePolynomials() async throws {
        let calculator = BodeCalculator()
        let entry = ExprEntry(text: "(s+1)^300")
        let result = try #require(await calculator.calculate(entries: [entry], start: -1, end: 1, points: 3))
        let middle = try #require(result[0].bodeData.first { $0.w == 1 })
        #expect(abs(middle.magDB - 300 * 10 * log10(2)) < 1e-8)

        let sum = ExprEntry(text: "s^300+1")
        let quotient = ExprEntry(text: "s^300/(s+2)", colorIndex: 1)
        let combined = try #require(await calculator.calculate(entries: [sum, quotient], start: -1, end: 1, points: 3))
        let sumAtOne = try #require(combined[0].bodeData.first { $0.w == 1 })
        let quotientAtOne = try #require(combined[1].bodeData.first { $0.w == 1 })
        #expect(abs(sumAtOne.magDB - 20 * log10(2)) < 1e-8)
        #expect(abs(quotientAtOne.magDB - (-10 * log10(5))) < 1e-8)

        // This number is larger than Int.max; entering it must not trap.
        let veryLarge = ExprEntry(text: "s^100000000000000000000")
        let extreme = try #require(await calculator.calculate(entries: [veryLarge], start: -1, end: 1, points: 3))
        #expect(extreme[0].errorMessage == nil)
        #expect(extreme[0].bodeData.contains { $0.w == 1 && abs($0.magDB) < 1e-8 })
    }

    @Test func allEightSeriesKeepTheirPointsWhenOneIsEdited() async throws {
        let calculator = BodeCalculator()
        var entries = (0..<8).map { ExprEntry(text: "\($0 + 1)/(s+1)", colorIndex: $0) }
        let before = try #require(await calculator.calculate(entries: entries, start: -2, end: 2, points: 3000))
        #expect(before.count == 8)
        #expect(before.allSatisfy { $0.bodeData.count == 3000 })

        entries[0].text = "2/(s+1)"
        let after = try #require(await calculator.calculate(entries: entries, start: -2, end: 2, points: 3000))
        #expect(after[0].bodeData[1500].magDB != before[0].bodeData[1500].magDB)
        for index in 1..<8 {
            #expect(after[index].id == before[index].id)
            #expect(after[index].bodeData.map(\.magDB) == before[index].bodeData.map(\.magDB))
            #expect(after[index].bodeData.map(\.phaseDeg) == before[index].bodeData.map(\.phaseDeg))
        }
    }

    @Test @MainActor func seriesAreDistinguishableWithoutColor() {
        let patterns = (0..<8).map {
            seriesStroke(index: $0, distinguishWithoutColor: true, width: 2.2).dash
        }
        #expect(Set(patterns).count == 8)
        #expect((0..<8).allSatisfy {
            seriesStroke(index: $0, distinguishWithoutColor: false, width: 2.2).dash.isEmpty
        })
    }

    @Test @MainActor func increasedContrastSeriesColorsRemainLegible() {
        func luminance(_ color: UIColor) -> Double {
            var red: CGFloat = 0, green: CGFloat = 0, blue: CGFloat = 0, alpha: CGFloat = 0
            #expect(color.getRed(&red, green: &green, blue: &blue, alpha: &alpha))
            func linear(_ value: CGFloat) -> Double {
                let value = Double(value)
                return value <= 0.04045 ? value / 12.92 : pow((value + 0.055) / 1.055, 2.4)
            }
            return 0.2126 * linear(red) + 0.7152 * linear(green) + 0.0722 * linear(blue)
        }
        func ratio(_ first: Double, _ second: Double) -> Double {
            (max(first, second) + 0.05) / (min(first, second) + 0.05)
        }

        let backgrounds = [UIUserInterfaceStyle.light, .dark].map { style in
            let traits = UITraitCollection { traits in
                traits.userInterfaceStyle = style
                traits.accessibilityContrast = .high
            }
            return luminance(UIColor.secondarySystemGroupedBackground.resolvedColor(with: traits))
        }
        for index in 0..<8 {
            let color = luminance(UIColor(seriesColor(index: index, highContrast: true)))
            #expect(ratio(color, 1) >= 4.5, "Badge text must contrast with its background")
            for background in backgrounds {
                #expect(ratio(color, background) >= 3, "Chart lines must contrast with their card")
            }
        }
    }

}
