//
//  BodePlotAppTests.swift
//  BodePlotAppTests
//
//  Created by Yuta Nisimatsu on 2025/12/23.
//

import Testing
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

}
