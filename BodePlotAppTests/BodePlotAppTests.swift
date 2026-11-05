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
