//
//  ContentView.swift
//  BodePlotApp
//
//  Created by Yuta Nisimatsu on 2025/12/23.
//

import SwiftUI
import Charts
import StoreKit

// =====================================================
//  Complex（Bode計算用）
// =====================================================
struct Complex: Equatable {
    var re: Double
    var im: Double
    static let zero = Complex(re: 0, im: 0)

    static func + (lhs: Complex, rhs: Complex) -> Complex {
        .init(re: lhs.re + rhs.re, im: lhs.im + rhs.im)
    }
    static func * (lhs: Complex, rhs: Complex) -> Complex {
        .init(re: lhs.re * rhs.re - lhs.im * rhs.im,
              im: lhs.re * rhs.im + lhs.im * rhs.re)
    }
    static func / (lhs: Complex, rhs: Complex) -> Complex {
        let denom = rhs.re * rhs.re + rhs.im * rhs.im
        if denom == 0 { return .init(re: .nan, im: .nan) }
        return .init(re: (lhs.re * rhs.re + lhs.im * rhs.im) / denom,
                     im: (lhs.im * rhs.re - lhs.re * rhs.im) / denom)
    }

    func abs() -> Double { hypot(re, im) }
    func phaseDeg() -> Double { atan2(im, re) * 180.0 / .pi }
}

extension Complex {
    func exp() -> Complex {
        let ea = Foundation.exp(re)
        return Complex(re: ea * cos(im), im: ea * sin(im))
    }

    func log() -> Complex {
        return Complex(re: Foundation.log(self.abs()), im: atan2(im, re))
    }

    func pow(_ p: Double) -> Complex {
        let l = self.log()
        return Complex(re: l.re * p, im: l.im * p).exp()
    }
}

// =====================================================
//  Bodeモデル
// =====================================================
struct BodePoint: Identifiable {
    var id: Double { w }
    let w: Double
    let logW: Double
    let magDB: Double
    let phaseDeg: Double
}

// =====================================================
//  Bode計算
// =====================================================
enum BodeCalc {
    static func logspace(startExp: Double, endExp: Double, count: Int) -> [Double] {
        guard count >= 2 else { return [pow(10, startExp)] }
        return (0..<count).map { i in
            let t = Double(i) / Double(count - 1)
            return pow(10, startExp + (endExp - startExp) * t)
        }
    }

    static func evalPoly(_ coeffsHighFirst: [Double], at s: Complex) -> Complex {
        var y = Complex.zero
        for c in coeffsHighFirst {
            y = y * s + Complex(re: c, im: 0)
        }
        return y
    }

    static func unwrapPhaseDeg(_ phases: [Double]) -> [Double] {
        guard !phases.isEmpty else { return [] }
        var out = [phases[0]]
        var offset = 0.0
        for i in 1..<phases.count {
            let prev = phases[i - 1] + offset
            var cur = phases[i] + offset
            let diff = cur - prev
            if diff > 180 { offset -= 360 }
            else if diff < -180 { offset += 360 }
            cur = phases[i] + offset
            out.append(cur)
        }
        return out
    }

    static func bode(numHighFirst: [Double], denHighFirst: [Double],
                     wStartExp: Double, wEndExp: Double, points: Int) -> [BodePoint] {

        let ws = logspace(startExp: wStartExp, endExp: wEndExp, count: points)

        var mags: [Double] = []
        mags.reserveCapacity(ws.count)

        var phases: [Double] = []
        phases.reserveCapacity(ws.count)

        for w in ws {
            let s = Complex(re: 0, im: w)
            let n = evalPoly(numHighFirst, at: s)
            let d = evalPoly(denHighFirst, at: s)
            let h = n / d

            mags.append(20.0 * log10(h.abs()))
            phases.append(h.phaseDeg())
        }

        let unwrapped = unwrapPhaseDeg(phases)

        return zip(zip(ws, mags), unwrapped).map { pair, ph in
            let (w, magDB) = pair
            return BodePoint(w: w, logW: log10(w), magDB: magDB, phaseDeg: ph)
        }
    }

    static func bodeEval(wStartExp: Double, wEndExp: Double, points: Int,
                         eval: (Complex) -> Complex) -> [BodePoint] {

        let ws = logspace(startExp: wStartExp, endExp: wEndExp, count: points)

        var mags: [Double] = []
        var phases: [Double] = []
        mags.reserveCapacity(ws.count)
        phases.reserveCapacity(ws.count)

        for w in ws {
            let s = Complex(re: 0, im: w)
            let h = eval(s)

            mags.append(20.0 * log10(h.abs()))
            phases.append(h.phaseDeg())
        }

        let unwrapped = unwrapPhaseDeg(phases)

        return zip(zip(ws, mags), unwrapped).map { pair, ph in
            let (w, magDB) = pair
            return BodePoint(w: w, logW: log10(w), magDB: magDB, phaseDeg: ph)
        }
    }
}

// =====================================================
//  ローカライズ用ヘルパー
// =====================================================
fileprivate func locError(_ key: String.LocalizationValue) -> ParseError {
    .message(String(localized: key))
}

// =====================================================
//  式 → (分子/分母) 多項式係数 パーサ
// =====================================================

typealias Poly = [Double]

fileprivate let polyEps = 1e-12

fileprivate func trimPoly(_ p: Poly) -> Poly {
    var p = p
    while p.count > 1, abs(p.last ?? 0) < polyEps { p.removeLast() }
    if p.isEmpty { return [0] }
    return p
}

fileprivate func polyAdd(_ a: Poly, _ b: Poly) -> Poly {
    let n = max(a.count, b.count)
    var out = Array(repeating: 0.0, count: n)
    for i in 0..<n {
        let av = i < a.count ? a[i] : 0
        let bv = i < b.count ? b[i] : 0
        out[i] = av + bv
    }
    return trimPoly(out)
}

fileprivate func polySub(_ a: Poly, _ b: Poly) -> Poly {
    let n = max(a.count, b.count)
    var out = Array(repeating: 0.0, count: n)
    for i in 0..<n {
        let av = i < a.count ? a[i] : 0
        let bv = i < b.count ? b[i] : 0
        out[i] = av - bv
    }
    return trimPoly(out)
}

fileprivate func polyMul(_ a: Poly, _ b: Poly) -> Poly {
    let a = trimPoly(a), b = trimPoly(b)
    if a.count == 1, abs(a[0]) < polyEps { return [0] }
    if b.count == 1, abs(b[0]) < polyEps { return [0] }
    var out = Array(repeating: 0.0, count: a.count + b.count - 1)
    for i in 0..<a.count {
        for j in 0..<b.count {
            out[i + j] += a[i] * b[j]
        }
    }
    return trimPoly(out)
}

fileprivate func polyScale(_ a: Poly, _ k: Double) -> Poly {
    trimPoly(a.map { $0 * k })
}

fileprivate func polyPow(_ a: Poly, _ e: Int) -> Poly {
    precondition(e >= 0)
    if e == 0 { return [1] }
    if e == 1 { return trimPoly(a) }
    var base = trimPoly(a)
    var exp = e
    var res: Poly = [1]
    while exp > 0 {
        if (exp & 1) == 1 { res = polyMul(res, base) }
        exp >>= 1
        if exp > 0 { base = polyMul(base, base) }
    }
    return trimPoly(res)
}

struct Rational {
    var num: Poly
    var den: Poly

    init(num: Poly, den: Poly) throws {
        let dn = trimPoly(den)
        if dn.count == 1, abs(dn[0]) < polyEps {
            throw locError("Denominator is zero")
        }
        self.num = trimPoly(num)
        self.den = dn
        normalize()
    }

    mutating func normalize() {
        let lead = den.last ?? 1
        if abs(lead) > polyEps {
            num = polyScale(num, 1.0 / lead)
            den = polyScale(den, 1.0 / lead)
        }
        if (den.last ?? 1) < 0 {
            num = polyScale(num, -1)
            den = polyScale(den, -1)
        }
    }

    static func + (lhs: Rational, rhs: Rational) throws -> Rational {
        let ad = polyMul(lhs.num, rhs.den)
        let cb = polyMul(rhs.num, lhs.den)
        let n = polyAdd(ad, cb)
        let d = polyMul(lhs.den, rhs.den)
        return try Rational(num: n, den: d)
    }
    static func - (lhs: Rational, rhs: Rational) throws -> Rational {
        let ad = polyMul(lhs.num, rhs.den)
        let cb = polyMul(rhs.num, lhs.den)
        let n = polySub(ad, cb)
        let d = polyMul(lhs.den, rhs.den)
        return try Rational(num: n, den: d)
    }
    static func * (lhs: Rational, rhs: Rational) throws -> Rational {
        return try Rational(num: polyMul(lhs.num, rhs.num),
                            den: polyMul(lhs.den, rhs.den))
    }
    static func / (lhs: Rational, rhs: Rational) throws -> Rational {
        return try Rational(num: polyMul(lhs.num, rhs.den),
                            den: polyMul(lhs.den, rhs.num))
    }

    func pow(_ e: Int) throws -> Rational {
        if e < 0 {
            throw locError("^ only supports non-negative integer exponents")
        }
        return try Rational(num: polyPow(num, e), den: polyPow(den, e))
    }

    func toHighFirstCoeffs() -> (num: [Double], den: [Double]) {
        func toHigh(_ p: Poly) -> [Double] {
            let p = trimPoly(p)
            var h = Array(p.reversed())
            while h.count > 1, abs(h.first ?? 0) < polyEps { h.removeFirst() }
            return h.isEmpty ? [0] : h
        }
        return (toHigh(num), toHigh(den))
    }
}

enum ParseError: Error, LocalizedError {
    case message(String)
    var errorDescription: String? {
        switch self { case .message(let s): return s }
    }
}

enum Token: Equatable {
    case number(Double)
    case ident(String)
    case op(Character)
    case lparen
    case rparen
    case end
}

struct Lexer {
    let input: String
    var i: String.Index

    init(_ raw: String) {
        let normalized = Lexer.normalize(raw)
        if let eq = normalized.firstIndex(of: "=") {
            self.input = String(normalized[normalized.index(after: eq)...])
        } else {
            self.input = normalized
        }
        self.i = self.input.startIndex
    }

    private static func normalize(_ raw: String) -> String {
        var out = ""
        out.reserveCapacity(raw.count)

        for ch in raw {
            if ch.isWhitespace { continue }

            switch ch {
            case "×", "∙", "·", "•", "・", "＊":
                out.append("*")
            case "÷", "／":
                out.append("/")
            case "−", "–", "—":
                out.append("-")
            case "＋":
                out.append("+")
            case "＾":
                out.append("^")
            case "（", "【", "［", "〔", "｛", "{", "[":
                out.append("(")
            case "）", "】", "］", "〕", "｝", "}", "]":
                out.append(")")
            default:
                out.append(ch)
            }
        }
        return out
    }

    mutating func nextToken() throws -> Token {
        while i < input.endIndex, input[i].isWhitespace { i = input.index(after: i) }
        if i >= input.endIndex { return .end }

        let ch = input[i]

        if ch == "(" { i = input.index(after: i); return .lparen }
        if ch == ")" { i = input.index(after: i); return .rparen }

        if "+-*/^".contains(ch) {
            i = input.index(after: i)
            return .op(ch)
        }

        if ch.isNumber || ch == "." {
            var j = i
            var sawDot = false
            while j < input.endIndex {
                let c = input[j]
                if c.isNumber { j = input.index(after: j); continue }
                if c == "." && !sawDot { sawDot = true; j = input.index(after: j); continue }
                break
            }
            let s = String(input[i..<j])
            i = j
            if let v = Double(s) {
                return .number(v)
            } else {
                throw locError("Cannot read number: \(s)")
            }
        }

        if ch.isLetter {
            var j = i
            while j < input.endIndex, input[j].isLetter { j = input.index(after: j) }
            let name = String(input[i..<j])
            i = j
            return .ident(name)
        }

        throw locError("Invalid character: \(String(ch))")
    }
}

final class Parser {
    private var tokens: [Token] = []
    private var pos: Int = 0

    init(_ expr: String) throws {
        var lx = Lexer(expr)
        var t: Token = try lx.nextToken()
        while t != .end {
            tokens.append(t)
            t = try lx.nextToken()
        }
        tokens.append(.end)
    }

    private func peek() -> Token { tokens[pos] }
    private func advance() { pos = min(pos + 1, tokens.count - 1) }

    private func matchOp(_ c: Character) -> Bool {
        if case .op(let x) = peek(), x == c {
            advance(); return true
        }
        return false
    }

    private func isImplicitMulStart(_ t: Token) -> Bool {
        switch t {
        case .number, .ident, .lparen: return true
        default: return false
        }
    }

    func parse() throws -> Rational {
        let r = try parseExpression()
        if peek() != .end {
            throw locError("Unexpected tokens after expression")
        }
        return r
    }

    private func parseExpression() throws -> Rational {
        var lhs = try parseTerm()
        while true {
            if matchOp("+") {
                let rhs = try parseTerm()
                lhs = try (lhs + rhs)
            } else if matchOp("-") {
                let rhs = try parseTerm()
                lhs = try (lhs - rhs)
            } else { break }
        }
        return lhs
    }

    private func parseTerm() throws -> Rational {
        var lhs = try parsePower()
        while true {
            if matchOp("*") {
                let rhs = try parsePower()
                lhs = try (lhs * rhs)
            } else if matchOp("/") {
                let rhs = try parsePower()
                lhs = try (lhs / rhs)
            } else if isImplicitMulStart(peek()) {
                let rhs = try parsePower()
                lhs = try (lhs * rhs)
            } else { break }
        }
        return lhs
    }

    private func parsePower() throws -> Rational {
        var base = try parseUnary()
        while matchOp("^") {
            let expTok = peek()
            guard case .number(let v) = expTok else {
                throw locError("^ must be followed by an integer (e.g. ^2)")
            }
            advance()
            let e = Int(v)
            if abs(Double(e) - v) > 1e-9 || e < 0 {
                throw locError("^ only supports non-negative integers (got: \(String(v)))")
            }
            base = try base.pow(e)
        }
        return base
    }

    private func parseUnary() throws -> Rational {
        if matchOp("+") { return try parseUnary() }
        if matchOp("-") {
            let r = try parseUnary()
            return try Rational(num: polyScale(r.num, -1), den: r.den)
        }
        return try parsePrimary()
    }

    private func parsePrimary() throws -> Rational {
        switch peek() {
        case .number(let v):
            advance()
            return try Rational(num: [v], den: [1])
        case .ident(let name):
            advance()
            if name.lowercased() != "s" {
                throw locError("Only variable 's' is supported (got: \(name))")
            }
            return try Rational(num: [0, 1], den: [1])
        case .lparen:
            advance()
            let r = try parseExpression()
            if peek() != .rparen {
                throw locError("Missing closing parenthesis ')'")
            }
            advance()
            return r
        case .end:
            throw locError("Please enter a valid expression")
        default:
            throw locError("Unparseable token: \(String(describing: peek()))")
        }
    }
}

// =====================================================
//  Complex式パーサ（非整数冪OK）
// =====================================================

indirect enum CExpr {
    case num(Double)
    case s
    case add(CExpr, CExpr)
    case sub(CExpr, CExpr)
    case mul(CExpr, CExpr)
    case div(CExpr, CExpr)
    case pow(CExpr, Double)
    case neg(CExpr)
    case sqrt(CExpr)

    func eval(_ sVal: Complex) -> Complex {
        switch self {
        case .num(let v): return Complex(re: v, im: 0)
        case .s: return sVal
        case .add(let a, let b): return a.eval(sVal) + b.eval(sVal)
        case .sub(let a, let b): return a.eval(sVal) + Complex(re: -b.eval(sVal).re, im: -b.eval(sVal).im)
        case .mul(let a, let b): return a.eval(sVal) * b.eval(sVal)
        case .div(let a, let b): return a.eval(sVal) / b.eval(sVal)
        case .neg(let x):
            let v = x.eval(sVal)
            return Complex(re: -v.re, im: -v.im)
        case .pow(let base, let p):
            return base.eval(sVal).pow(p)
        case .sqrt(let x):
            return x.eval(sVal).pow(0.5)
        }
    }
}

final class ComplexExprParser {
    private var tokens: [Token] = []
    private var pos: Int = 0

    init(_ expr: String) throws {
        var lx = Lexer(expr)
        var t: Token = try lx.nextToken()
        while t != .end {
            tokens.append(t)
            t = try lx.nextToken()
        }
        tokens.append(.end)
    }

    private func peek() -> Token { tokens[pos] }
    private func advance() { pos = min(pos + 1, tokens.count - 1) }

    private func matchOp(_ c: Character) -> Bool {
        if case .op(let x) = peek(), x == c { advance(); return true }
        return false
    }

    private func isImplicitMulStart(_ t: Token) -> Bool {
        switch t {
        case .number, .ident, .lparen: return true
        default: return false
        }
    }

    func parse() throws -> CExpr {
        let e = try parseExpression()
        if peek() != .end {
            throw locError("Unexpected tokens after expression")
        }
        return e
    }

    private func parseExpression() throws -> CExpr {
        var lhs = try parseTerm()
        while true {
            if matchOp("+") {
                let rhs = try parseTerm()
                lhs = .add(lhs, rhs)
            } else if matchOp("-") {
                let rhs = try parseTerm()
                lhs = .sub(lhs, rhs)
            } else { break }
        }
        return lhs
    }

    private func parseTerm() throws -> CExpr {
        var lhs = try parsePower()
        while true {
            if matchOp("*") {
                let rhs = try parsePower()
                lhs = .mul(lhs, rhs)
            } else if matchOp("/") {
                let rhs = try parsePower()
                lhs = .div(lhs, rhs)
            } else if isImplicitMulStart(peek()) {
                let rhs = try parsePower()
                lhs = .mul(lhs, rhs)
            } else { break }
        }
        return lhs
    }

    private func parsePower() throws -> CExpr {
        var base = try parseUnary()
        while matchOp("^") {
            let p = try parseRealConst()
            base = .pow(base, p)
        }
        return base
    }

    private func parseUnary() throws -> CExpr {
        if matchOp("+") { return try parseUnary() }
        if matchOp("-") { return .neg(try parseUnary()) }
        return try parsePrimary()
    }

    private func parsePrimary() throws -> CExpr {
        switch peek() {
        case .number(let v):
            advance()
            return .num(v)

        case .ident(let name):
            advance()
            let lower = name.lowercased()
            if lower == "s" { return .s }
            if lower == "sqrt" {
                guard peek() == .lparen else {
                    throw locError("sqrt must be followed by ()")
                }
                advance()
                let arg = try parseExpression()
                if peek() != .rparen {
                    throw locError("Missing closing parenthesis ')'")
                }
                advance()
                return .sqrt(arg)
            }
            throw locError("Only variable 's' and function sqrt() are supported (got: \(name))")

        case .lparen:
            advance()
            let e = try parseExpression()
            if peek() != .rparen {
                throw locError("Missing closing parenthesis ')'")
            }
            advance()
            return e

        case .end:
            throw locError("Please enter a valid expression")
        default:
            throw locError("Unparseable token: \(String(describing: peek()))")
        }
    }

    private func parseRealConst() throws -> Double {
        if peek() == .lparen {
            advance()
            let v = try parseRealExpr()
            if peek() != .rparen {
                throw locError("Missing closing parenthesis ')' in exponent")
            }
            advance()
            return v
        }
        return try parseRealExpr()
    }

    private func parseRealExpr() throws -> Double {
        var lhs = try parseRealTerm()
        while true {
            if matchOp("+") { lhs += try parseRealTerm() }
            else if matchOp("-") { lhs -= try parseRealTerm() }
            else { break }
        }
        return lhs
    }

    private func parseRealTerm() throws -> Double {
        var lhs = try parseRealUnary()
        while true {
            if matchOp("*") { lhs *= try parseRealUnary() }
            else if matchOp("/") { lhs /= try parseRealUnary() }
            else { break }
        }
        return lhs
    }

    private func parseRealUnary() throws -> Double {
        if matchOp("+") { return try parseRealUnary() }
        if matchOp("-") { return -(try parseRealUnary()) }
        return try parseRealPrimary()
    }

    private func parseRealPrimary() throws -> Double {
        switch peek() {
        case .number(let v):
            advance(); return v
        case .lparen:
            advance()
            let v = try parseRealExpr()
            if peek() != .rparen {
                throw locError("Missing closing parenthesis ')' in exponent")
            }
            advance()
            return v
        case .ident:
            throw locError("Cannot use 's' in exponent (must be a real constant)")
        default:
            throw locError("Cannot read exponent")
        }
    }
}

// =====================================================
//  レビュー促進
// =====================================================
struct ReviewPromptModifier: ViewModifier {
    @AppStorage("bode_firstLaunchTS")     private var firstLaunchTS: Double = 0
    @AppStorage("bode_reviewPromptDone")  private var promptDone: Bool = false

    @State private var showAlert = false
    @Environment(\.requestReview) private var requestReview

    private let oneWeek: TimeInterval = 7 * 24 * 60 * 60

    func body(content: Content) -> some View {
        content
            .task { await checkPrompt() }
            .alert("Enjoying Bode Plot?", isPresented: $showAlert) {
                Button("Write a Review") {
                    promptDone = true
                    requestReview()
                }
                Button("Maybe Later") {
                    firstLaunchTS = Date().timeIntervalSince1970
                }
                Button("No Thanks", role: .cancel) {
                    promptDone = true
                }
            } message: {
                Text("If you find this app useful, please consider leaving a review on the App Store. It really helps!")
            }
    }

    private func checkPrompt() async {
        if firstLaunchTS == 0 {
            firstLaunchTS = Date().timeIntervalSince1970
            return
        }
        if promptDone { return }
        let elapsed = Date().timeIntervalSince1970 - firstLaunchTS
        guard elapsed >= oneWeek else { return }
        try? await Task.sleep(nanoseconds: 1_500_000_000)
        await MainActor.run { showAlert = true }
    }
}

extension View {
    func reviewPrompt() -> some View { modifier(ReviewPromptModifier()) }
}

// =====================================================
//  モデル
// =====================================================
struct ExprEntry: Codable, Identifiable, Equatable {
    let id: UUID
    var text: String
    var colorIndex: Int

    init(id: UUID = UUID(), text: String = "", colorIndex: Int = 0) {
        self.id = id
        self.text = text
        self.colorIndex = colorIndex
    }
}

struct BodeSeries: Identifiable {
    let id: UUID
    let color: Color
    let label: String
    let bodeData: [BodePoint]
    let errorMessage: String?
}

private let seriesColors: [Color] = [
    Color(red: 0.00, green: 0.48, blue: 1.00),
    Color(red: 1.00, green: 0.58, blue: 0.00),
    Color(red: 0.20, green: 0.78, blue: 0.35),
    Color(red: 1.00, green: 0.23, blue: 0.19),
    Color(red: 0.69, green: 0.32, blue: 0.87),
    Color(red: 0.00, green: 0.78, blue: 0.85),
    Color(red: 1.00, green: 0.18, blue: 0.59),
    Color(red: 0.35, green: 0.34, blue: 0.84),
]

private func subscriptDigits(_ n: Int) -> String {
    let map: [Character: Character] = [
        "0":"₀","1":"₁","2":"₂","3":"₃","4":"₄",
        "5":"₅","6":"₆","7":"₇","8":"₈","9":"₉"
    ]
    return String(String(n).map { map[$0] ?? $0 })
}

// 数式プレースホルダや軸ラベルは数学記号のため翻訳不要
private let exprPlaceholder = "G(s) = 100 / (s(s+1)^2)"
private let gainAxisLabel   = "20 log |G(jω)|  [dB]"
private let phaseAxisLabel  = "∠G(jω)  [deg]"

// =====================================================
//  小さなビュー部品
// =====================================================
struct SeriesBadge: View {
    let index: Int
    let color: Color
    var size: CGFloat = 28

    var body: some View {
        ZStack {
            Circle()
                .fill(color)
            Text(verbatim: "\(index)")
                .font(.system(size: size * 0.45, weight: .bold, design: .rounded))
                .foregroundStyle(.white)
        }
        .frame(width: size, height: size)
        .shadow(color: color.opacity(0.3), radius: 2, y: 1)
    }
}

struct ChartLegendView: View {
    struct Item: Identifiable {
        let id: UUID
        let label: String
        let color: Color
    }
    let items: [Item]

    var body: some View {
        ScrollView(.horizontal, showsIndicators: false) {
            HStack(spacing: 14) {
                ForEach(items) { item in
                    HStack(spacing: 6) {
                        RoundedRectangle(cornerRadius: 1.5)
                            .fill(item.color)
                            .frame(width: 18, height: 3)
                        Text(verbatim: item.label)
                            .font(.caption.weight(.medium))
                            .foregroundStyle(.primary)
                    }
                }
            }
            .padding(.horizontal, 2)
        }
    }
}

struct EmptyChartHintView: View {
    var body: some View {
        VStack(spacing: 12) {
            Image(systemName: "waveform.path")
                .font(.system(size: 38, weight: .light))
                .foregroundStyle(.tertiary)
            Text("Enter a transfer function above to plot its Bode diagram here")
                .font(.callout)
                .foregroundStyle(.secondary)
                .multilineTextAlignment(.center)
                .padding(.horizontal, 24)
        }
        .frame(maxWidth: .infinity)
        .frame(height: 200)
    }
}

// =====================================================
//  1行分の式入力ビュー
// =====================================================
struct ExprRowView: View {
    let index: Int
    @Binding var text: String
    let color: Color
    let error: String?
    let canDelete: Bool
    let onDelete: () -> Void
    @FocusState.Binding var focusedID: UUID?
    let entryID: UUID

    var body: some View {
        VStack(alignment: .leading, spacing: 6) {
            HStack(alignment: .top, spacing: 12) {
                SeriesBadge(index: index, color: color)
                    .padding(.top, 4)

                TextField(exprPlaceholder, text: $text, axis: .vertical)
                    .textFieldStyle(.plain)
                    .font(.body.monospaced())
                    .padding(.horizontal, 12)
                    .padding(.vertical, 10)
                    .background(
                        RoundedRectangle(cornerRadius: 10, style: .continuous)
                            .fill(Color(.tertiarySystemBackground))
                    )
                    .overlay(
                        RoundedRectangle(cornerRadius: 10, style: .continuous)
                            .strokeBorder(
                                focusedID == entryID ? color.opacity(0.6)
                                                     : Color.secondary.opacity(0.15),
                                lineWidth: focusedID == entryID ? 1.5 : 1
                            )
                    )
                    .lineLimit(1...4)
                    .focused($focusedID, equals: entryID)
                    .submitLabel(.done)
                    .autocorrectionDisabled(true)
                    .textInputAutocapitalization(.never)

                if canDelete {
                    Button(role: .destructive, action: onDelete) {
                        Image(systemName: "minus.circle.fill")
                            .font(.title3)
                            .foregroundStyle(.red.opacity(0.85))
                            .symbolRenderingMode(.hierarchical)
                    }
                    .buttonStyle(.plain)
                    .padding(.top, 6)
                    .accessibilityLabel("Delete")
                }
            }

            if let error {
                HStack(alignment: .top, spacing: 6) {
                    Image(systemName: "exclamationmark.triangle.fill")
                        .font(.caption2)
                        .symbolRenderingMode(.hierarchical)
                    Text(verbatim: error)
                        .font(.caption)
                }
                .foregroundStyle(.red)
                .padding(.leading, 40)
            }
        }
    }
}

// =====================================================
//  メインビュー
// =====================================================
struct ContentView: View {

    @AppStorage("bode_exprsJSON")      private var exprsJSON: String = ""
    @AppStorage("bode_wStartExp")      private var wStartExp: Double = -2
    @AppStorage("bode_wEndExp")        private var wEndExp: Double = 2
    @AppStorage("bode_points")         private var points: Double = 1000
    @AppStorage("bode_isAdvancedOpen") private var isAdvancedOpen: Bool = false

    @State private var entries: [ExprEntry] = [ExprEntry()]
    @State private var allSeries: [BodeSeries] = []
    @State private var recalcTask: Task<Void, Never>?
    @FocusState private var focusedID: UUID?

    private var activeSeries: [BodeSeries] {
        allSeries.filter { !$0.bodeData.isEmpty }
    }

    private var canAddMore: Bool { entries.count < seriesColors.count }

    var body: some View {
        NavigationStack {
            ScrollView {
                LazyVStack(alignment: .leading, spacing: 20) {

                    expressionSection
                    chartSection
                    advancedSection

                    Text("Horizontal axis is angular frequency ω on a logarithmic scale (log₁₀ω).")
                        .font(.footnote)
                        .foregroundStyle(.secondary)
                        .padding(.top, 4)
                }
                .padding()
            }
            .background(Color(.systemGroupedBackground).ignoresSafeArea())
            .navigationTitle("Bode Plot")
            .navigationBarTitleDisplayMode(.inline)
            .scrollDismissesKeyboard(.interactively)
            .toolbar {
                ToolbarItemGroup(placement: .keyboard) {
                    Spacer()
                    Button {
                        focusedID = nil
                    } label: {
                        Text("Done").bold()
                    }
                }
            }
        }
        .onAppear {
            loadEntries()
            scheduleRecalc()
        }
        .onChange(of: entries) { _, _ in
            saveEntries()
            scheduleRecalc()
        }
        .onChange(of: wStartExp) { _, _ in scheduleRecalc() }
        .onChange(of: wEndExp)   { _, _ in scheduleRecalc() }
        .onChange(of: points)    { _, _ in scheduleRecalc() }
        .reviewPrompt()
    }

    // MARK: - Sections

    @ViewBuilder
    private var expressionSection: some View {
        VStack(alignment: .leading, spacing: 0) {
            SectionHeader(titleKey: "Transfer Functions", systemImage: "function")

            VStack(spacing: 0) {
                ForEach(Array(entries.enumerated()), id: \.element.id) { idx, entry in
                    ExprRowView(
                        index: idx + 1,
                        text: $entries[idx].text,
                        color: seriesColors[entry.colorIndex % seriesColors.count],
                        error: allSeries.first(where: { $0.id == entry.id })?.errorMessage,
                        canDelete: entries.count > 1,
                        onDelete: { deleteEntry(at: idx) },
                        focusedID: $focusedID,
                        entryID: entry.id
                    )
                    .padding(.horizontal, 14)
                    .padding(.vertical, 10)
                    .transition(.asymmetric(
                        insertion: .opacity.combined(with: .move(edge: .top)),
                        removal: .opacity.combined(with: .scale(scale: 0.95))
                    ))

                    if idx < entries.count - 1 {
                        Divider().padding(.leading, 14)
                    }
                }

                Divider()

                HStack(spacing: 12) {
                    Button {
                        withAnimation(.spring(response: 0.35, dampingFraction: 0.85)) {
                            addEntry()
                        }
                    } label: {
                        Label("Add Expression", systemImage: "plus.circle.fill")
                            .font(.callout.weight(.semibold))
                    }
                    .buttonStyle(.borderless)
                    .disabled(!canAddMore)

                    Spacer()

                    Button(role: .destructive) {
                        withAnimation(.easeInOut(duration: 0.25)) {
                            entries = [ExprEntry(colorIndex: 0)]
                            focusedID = nil
                        }
                    } label: {
                        Label("Clear All", systemImage: "trash")
                            .font(.callout)
                    }
                    .buttonStyle(.borderless)
                    .disabled(entries.count == 1 && (entries.first?.text.isEmpty ?? true))
                }
                .padding(.horizontal, 14)
                .padding(.vertical, 10)
            }
            .background(
                RoundedRectangle(cornerRadius: 14, style: .continuous)
                    .fill(Color(.secondarySystemGroupedBackground))
            )
        }
    }

    @ViewBuilder
    private var chartSection: some View {
        VStack(alignment: .leading, spacing: 14) {
            if activeSeries.isEmpty {
                EmptyChartHintView()
                    .background(
                        RoundedRectangle(cornerRadius: 14, style: .continuous)
                            .fill(Color(.secondarySystemGroupedBackground))
                    )
            } else {
                if activeSeries.count > 1 {
                    ChartLegendView(items: activeSeries.enumerated().map { idx, s in
                        ChartLegendView.Item(
                            id: s.id,
                            label: "G" + subscriptDigits(idx + 1),
                            color: s.color
                        )
                    })
                    .padding(.horizontal, 14)
                    .padding(.vertical, 10)
                    .background(
                        RoundedRectangle(cornerRadius: 12, style: .continuous)
                            .fill(Color(.secondarySystemGroupedBackground))
                    )
                }

                BodeChartsView(
                    series: activeSeries,
                    wStartExp: min(wStartExp, wEndExp - 0.5),
                    wEndExp: max(wEndExp, wStartExp + 0.5)
                )
            }
        }
    }

    @ViewBuilder
    private var advancedSection: some View {
        VStack(alignment: .leading, spacing: 0) {
            Button {
                withAnimation(.spring(response: 0.35, dampingFraction: 0.85)) {
                    isAdvancedOpen.toggle()
                }
            } label: {
                HStack(spacing: 8) {
                    Image(systemName: "slider.horizontal.3")
                        .font(.subheadline.weight(.semibold))
                        .foregroundStyle(.secondary)
                    Text("Advanced Settings")
                        .font(.subheadline.weight(.semibold))
                        .foregroundStyle(.primary)
                    Spacer()
                    Image(systemName: "chevron.right")
                        .font(.caption.weight(.semibold))
                        .foregroundStyle(.secondary)
                        .rotationEffect(.degrees(isAdvancedOpen ? 90 : 0))
                }
                .padding(.horizontal, 14)
                .padding(.vertical, 14)
                .contentShape(Rectangle())
            }
            .buttonStyle(.plain)

            if isAdvancedOpen {
                Divider()
                VStack(alignment: .leading, spacing: 16) {
                    Text("Frequency Range (rad/s)")
                        .font(.caption.weight(.semibold))
                        .foregroundStyle(.secondary)
                        .padding(.top, 4)

                    sliderRow(labelKey: "Start 10^",
                              value: $wStartExp,
                              range: -6...3, step: 0.5,
                              format: "%.1f", width: 50)

                    sliderRow(labelKey: "End 10^",
                              value: $wEndExp,
                              range: -4...6, step: 0.5,
                              format: "%.1f", width: 50)

                    sliderRow(labelKey: "Points",
                              value: $points,
                              range: 200...3000, step: 100,
                              format: "%.0f", width: 60)

                    Button {
                        withAnimation { resetAdvancedSettings() }
                    } label: {
                        Label("Reset to Defaults", systemImage: "arrow.counterclockwise")
                            .font(.callout.weight(.medium))
                            .frame(maxWidth: .infinity)
                    }
                    .buttonStyle(.bordered)
                    .controlSize(.regular)
                    .padding(.top, 6)
                }
                .padding(.horizontal, 14)
                .padding(.bottom, 14)
                .transition(.opacity.combined(with: .move(edge: .top)))
            }
        }
        .background(
            RoundedRectangle(cornerRadius: 14, style: .continuous)
                .fill(Color(.secondarySystemGroupedBackground))
        )
    }

    @ViewBuilder
    private func sliderRow(labelKey: LocalizedStringKey,
                           value: Binding<Double>,
                           range: ClosedRange<Double>,
                           step: Double,
                           format: String,
                           width: CGFloat) -> some View {
        HStack(spacing: 12) {
            Text(labelKey)
                .font(.subheadline)
                .frame(width: 78, alignment: .leading)
            Slider(value: value, in: range, step: step)
            Text(verbatim: String(format: format, value.wrappedValue))
                .font(.subheadline.monospacedDigit())
                .foregroundStyle(.secondary)
                .frame(width: width, alignment: .trailing)
        }
    }

    // MARK: - Entry management

    private func addEntry() {
        let used = Set(entries.map { $0.colorIndex })
        let next = (0..<seriesColors.count).first { !used.contains($0) }
            ?? (entries.count % seriesColors.count)
        entries.append(ExprEntry(colorIndex: next))
    }

    private func deleteEntry(at idx: Int) {
        guard entries.indices.contains(idx) else { return }
        withAnimation(.spring(response: 0.35, dampingFraction: 0.85)) {
            entries.remove(at: idx)
            if entries.isEmpty { entries = [ExprEntry(colorIndex: 0)] }
        }
    }

    // MARK: - Persistence

    private func loadEntries() {
        if !exprsJSON.isEmpty,
           let data = exprsJSON.data(using: .utf8),
           let arr = try? JSONDecoder().decode([ExprEntry].self, from: data),
           !arr.isEmpty {
            entries = arr
            return
        }
        if let old = UserDefaults.standard.string(forKey: "bode_exprText"), !old.isEmpty {
            entries = [ExprEntry(text: old, colorIndex: 0)]
            saveEntries()
            return
        }
        entries = [ExprEntry(colorIndex: 0)]
    }

    private func saveEntries() {
        guard let data = try? JSONEncoder().encode(entries),
              let str = String(data: data, encoding: .utf8) else { return }
        exprsJSON = str
    }

    // MARK: - Recalc

    private func scheduleRecalc() {
        recalcTask?.cancel()
        recalcTask = Task {
            try? await Task.sleep(nanoseconds: 150_000_000)
            if Task.isCancelled { return }
            await MainActor.run { recalc() }
        }
    }

    private func recalc() {
        let start = min(wStartExp, wEndExp - 0.5)
        let end   = max(wEndExp,   wStartExp + 0.5)

        allSeries = entries.enumerated().map { idx, entry in
            let color = seriesColors[entry.colorIndex % seriesColors.count]
            let label = "G" + subscriptDigits(idx + 1)
            let trimmed = entry.text.trimmingCharacters(in: .whitespacesAndNewlines)

            guard !trimmed.isEmpty else {
                return BodeSeries(id: entry.id, color: color, label: label,
                                  bodeData: [], errorMessage: nil)
            }

            do {
                let parser = try Parser(entry.text)
                let r = try parser.parse()
                let coeffs = r.toHighFirstCoeffs()
                let raw = BodeCalc.bode(
                    numHighFirst: coeffs.num, denHighFirst: coeffs.den,
                    wStartExp: start, wEndExp: end, points: Int(points)
                )
                let filtered = raw.filter {
                    $0.logW.isFinite && $0.magDB.isFinite && $0.phaseDeg.isFinite
                }
                if filtered.isEmpty {
                    return BodeSeries(id: entry.id, color: color, label: label,
                                      bodeData: [],
                                      errorMessage: String(localized: "Result is NaN or infinity"))
                }
                return BodeSeries(id: entry.id, color: color, label: label,
                                  bodeData: filtered, errorMessage: nil)
            } catch {
                do {
                    let p2 = try ComplexExprParser(entry.text)
                    let expr = try p2.parse()
                    let raw = BodeCalc.bodeEval(
                        wStartExp: start, wEndExp: end, points: Int(points)
                    ) { s in expr.eval(s) }
                    let filtered = raw.filter {
                        $0.logW.isFinite && $0.magDB.isFinite && $0.phaseDeg.isFinite
                    }
                    if filtered.isEmpty {
                        return BodeSeries(id: entry.id, color: color, label: label,
                                          bodeData: [],
                                          errorMessage: String(localized: "Result is NaN or infinity"))
                    }
                    return BodeSeries(id: entry.id, color: color, label: label,
                                      bodeData: filtered, errorMessage: nil)
                } catch {
                    let msg = (error as? LocalizedError)?.errorDescription ?? "\(error)"
                    return BodeSeries(id: entry.id, color: color, label: label,
                                      bodeData: [], errorMessage: msg)
                }
            }
        }
    }

    private func resetAdvancedSettings() {
        wStartExp = -2
        wEndExp = 2
        points = 1000
    }
}

// =====================================================
//  セクションヘッダー
// =====================================================
struct SectionHeader: View {
    let titleKey: LocalizedStringKey
    let systemImage: String

    var body: some View {
        HStack(spacing: 8) {
            Image(systemName: systemImage)
                .font(.subheadline.weight(.semibold))
                .foregroundStyle(.secondary)
            Text(titleKey)
                .font(.subheadline.weight(.semibold))
                .textCase(.uppercase)
                .foregroundStyle(.secondary)
        }
        .padding(.horizontal, 6)
        .padding(.bottom, 6)
    }
}

// =====================================================
//  Charts（ゲイン/位相）
// =====================================================
struct BodeChartsView: View {
    let series: [BodeSeries]
    let wStartExp: Double
    let wEndExp: Double

    private var phaseMinTick: Double {
        let minV = series.flatMap { $0.bodeData }.map { $0.phaseDeg }.min() ?? -180
        return floor(minV / 90.0) * 90.0
    }

    private var phaseMaxTick: Double {
        let maxV = series.flatMap { $0.bodeData }.map { $0.phaseDeg }.max() ?? 180
        return ceil(maxV / 90.0) * 90.0
    }

    var body: some View {
        VStack(alignment: .leading, spacing: 16) {
            chartCard(titleKey: "Magnitude", subtitle: gainAxisLabel) {
                Chart {
                    ForEach(series) { s in
                        ForEach(s.bodeData.sorted { $0.logW < $1.logW }) { p in
                            LineMark(
                                x: .value("logW", p.logW),
                                y: .value("Mag(dB)", p.magDB),
                                series: .value("Series", s.id.uuidString)
                            )
                            .foregroundStyle(s.color)
                            .lineStyle(StrokeStyle(lineWidth: 2.2, lineCap: .round, lineJoin: .round))
                            .interpolationMethod(.catmullRom)
                        }
                    }
                }
                .frame(height: 240)
                .chartXScale(domain: wStartExp...wEndExp)
                .chartXAxis { bodeXAxis(startExp: wStartExp, endExp: wEndExp) }
                .chartYAxis { AxisMarks(position: .leading) }
                .chartLegend(.hidden)
            }

            chartCard(titleKey: "Phase", subtitle: phaseAxisLabel) {
                Chart {
                    ForEach(series) { s in
                        ForEach(s.bodeData.sorted { $0.logW < $1.logW }) { p in
                            LineMark(
                                x: .value("logW", p.logW),
                                y: .value("Phase(deg)", p.phaseDeg),
                                series: .value("Series", s.id.uuidString)
                            )
                            .foregroundStyle(s.color)
                            .lineStyle(StrokeStyle(lineWidth: 2.2, lineCap: .round, lineJoin: .round))
                            .interpolationMethod(.catmullRom)
                        }
                    }
                }
                .frame(height: 240)
                .chartXScale(domain: wStartExp...wEndExp)
                .chartXAxis { bodeXAxis(startExp: wStartExp, endExp: wEndExp) }
                .chartYAxis {
                    AxisMarks(
                        position: .leading,
                        values: stride(from: phaseMinTick, through: phaseMaxTick, by: 90.0).map { $0 }
                    ) { value in
                        AxisGridLine()
                        AxisTick()
                        AxisValueLabel {
                            if let y = value.as(Double.self) {
                                Text(verbatim: String(format: "%.0f°", y))
                            }
                        }
                    }
                }
                .chartLegend(.hidden)
            }
        }
    }

    @ViewBuilder
    private func chartCard<Content: View>(titleKey: LocalizedStringKey,
                                          subtitle: String,
                                          @ViewBuilder content: () -> Content) -> some View {
        VStack(alignment: .leading, spacing: 6) {
            HStack(alignment: .firstTextBaseline) {
                Text(titleKey)
                    .font(.headline)
                Spacer()
                Text(verbatim: subtitle)
                    .font(.caption2.monospaced())
                    .foregroundStyle(.secondary)
            }
            content()
        }
        .padding(14)
        .background(
            RoundedRectangle(cornerRadius: 14, style: .continuous)
                .fill(Color(.secondarySystemGroupedBackground))
        )
    }

    @AxisContentBuilder
    private func bodeXAxis(startExp: Double, endExp: Double) -> some AxisContent {
        let step = 1.0
        let exps: [Double] = {
            let start = ceil(startExp)
            let end = floor(endExp)
            if start > end {
                return [round((startExp + endExp) * 0.5)]
            } else {
                return stride(from: start, through: end, by: step).map { $0 }
            }
        }()

        AxisMarks(values: exps) { value in
            AxisGridLine()
            AxisTick()
            AxisValueLabel {
                if let x = value.as(Double.self) {
                    Text(verbatim: format10PowLabel(exp: x))
                }
            }
        }
    }

    private func superscriptInt(_ n: Int) -> String {
        let map: [Character: String] = [
            "0":"⁰","1":"¹","2":"²","3":"³","4":"⁴","5":"⁵","6":"⁶","7":"⁷","8":"⁸","9":"⁹",
            "-":"⁻"
        ]
        return String(String(n).map { map[$0] ?? String($0) }.joined())
    }

    private func format10PowLabel(exp: Double) -> String {
        let e = Int(round(exp))
        return "10" + superscriptInt(e)
    }
}
