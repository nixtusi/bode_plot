//
//  ContentView.swift
//  BodePlotApp
//
//  Created by Yuta Nisimatsu on 2025/12/23.
//

import SwiftUI
import Charts

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
        // exp(a+jb)=exp(a)(cos b + j sin b)
        let ea = Foundation.exp(re)
        return Complex(re: ea * cos(im), im: ea * sin(im))
    }

    func log() -> Complex {
        // principal value log(z)=ln|z| + j arg(z)
        return Complex(re: Foundation.log(self.abs()), im: atan2(im, re))
    }

    func pow(_ p: Double) -> Complex {
        // z^p = exp(p * log(z))  (principal branch)
        let l = self.log()
        return Complex(re: l.re * p, im: l.im * p).exp()
    }
}

// =====================================================
//  Bodeモデル
// =====================================================
struct BodePoint: Identifiable {
    var id: Double { w }   // ← w をIDに（並び順・描画安定化）
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

    /// Horner: coeffs are highest order first
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
            let s = Complex(re: 0, im: w) // jω
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
            let s = Complex(re: 0, im: w) // jω
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
//  ここから「式 → (分子/分母)多項式係数」パーサ
// =====================================================

typealias Poly = [Double] // 低次→高次（a0 + a1*s + ...）

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
            throw ParseError.message("分母が 0 になってる")
        }
        self.num = trimPoly(num)
        self.den = dn
        normalize()
    }

    mutating func normalize() {
        // 分母の最高次係数で軽く正規化（値が大きくなりすぎるの防止）
        let lead = den.last ?? 1
        if abs(lead) > polyEps {
            num = polyScale(num, 1.0 / lead)
            den = polyScale(den, 1.0 / lead)
        }
        // 分母の符号を揃える（最高次係数を正に）
        if (den.last ?? 1) < 0 {
            num = polyScale(num, -1)
            den = polyScale(den, -1)
        }
    }

    static func + (lhs: Rational, rhs: Rational) throws -> Rational {
        // a/b + c/d = (ad + cb) / bd
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
        // a/b ÷ c/d = ad/bc
        return try Rational(num: polyMul(lhs.num, rhs.den),
                            den: polyMul(lhs.den, rhs.num))
    }

    func pow(_ e: Int) throws -> Rational {
        if e < 0 { throw ParseError.message("^は0以上の整数だけ対応") }
        return try Rational(num: polyPow(num, e), den: polyPow(den, e))
    }

    func toHighFirstCoeffs() -> (num: [Double], den: [Double]) {
        // low-first -> high-first
        func toHigh(_ p: Poly) -> [Double] {
            let p = trimPoly(p)
            var h = Array(p.reversed())
            // 先頭の超微小を削る（高次側）
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
    case ident(String)     // "s"
    case op(Character)     // + - * / ^
    case lparen
    case rparen
    case end
}

struct Lexer {
    let input: String
    var i: String.Index

    init(_ raw: String) {
        // 入力ゆらぎに強くする（全角/記号ゆらぎ/波括弧など）
        // 例:
        //  - 10×9+s  -> 10*9+s
        //  - 1/{s(s+1)} -> 1/(s(s+1))
        //  - 全角括弧（ ）や｛ ｝もOK
        let normalized = Lexer.normalize(raw)

        // "G(s)=..." が来る前提：最初の '=' の右側だけ使う
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
            // すべての空白（全角スペース含む）を除去
            if ch.isWhitespace { continue }

            switch ch {
            // 掛け算
            case "×", "∙", "·", "•", "・", "＊":
                out.append("*")
            // 割り算
            case "÷", "／":
                out.append("/")
            // マイナス（Unicodeのダッシュ類）
            case "−", "–", "—":
                out.append("-")
            // プラス（全角）
            case "＋":
                out.append("+")
            // 冪（全角）
            case "＾":
                out.append("^")

            // 括弧：丸括弧以外も許容
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
                throw ParseError.message("数値が読めない: \(s)")
            }
        }

        if ch.isLetter {
            var j = i
            while j < input.endIndex, input[j].isLetter { j = input.index(after: j) }
            let name = String(input[i..<j])
            i = j
            return .ident(name)
        }

        throw ParseError.message("不正な文字: \(ch)")
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
        case .number, .ident, .lparen:
            return true
        default:
            return false
        }
    }

    func parse() throws -> Rational {
        let r = try parseExpression()
        if peek() != .end {
            throw ParseError.message("式の後ろに余計なものがある")
        }
        return r
    }

    // expression := term { (+|-) term }
    private func parseExpression() throws -> Rational {
        var lhs = try parseTerm()
        while true {
            if matchOp("+") {
                let rhs = try parseTerm()
                lhs = try (lhs + rhs)
            } else if matchOp("-") {
                let rhs = try parseTerm()
                lhs = try (lhs - rhs)
            } else {
                break
            }
        }
        return lhs
    }

    // term := power { (*|/|implicitMul) power }
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
                // 暗黙の掛け算: 2(s+1), (s+1)(s+2), s(s+1) など
                let rhs = try parsePower()
                lhs = try (lhs * rhs)
            } else {
                break
            }
        }
        return lhs
    }

    // power := unary { ^ INT }   （^ は 0以上整数のみ）
    private func parsePower() throws -> Rational {
        var base = try parseUnary()
        while matchOp("^") {
            let expTok = peek()
            guard case .number(let v) = expTok else {
                throw ParseError.message("^ の後ろは整数（例: ^2）だけ対応")
            }
            advance()
            let e = Int(v)
            if abs(Double(e) - v) > 1e-9 || e < 0 {
                throw ParseError.message("^ は 0以上の整数だけ対応（今: \(v)）")
            }
            base = try base.pow(e)
        }
        return base
    }

    // unary := (+ unary) | (- unary) | primary
    private func parseUnary() throws -> Rational {
        if matchOp("+") { return try parseUnary() }
        if matchOp("-") {
            let r = try parseUnary()
            return try Rational(num: polyScale(r.num, -1), den: r.den)
        }
        return try parsePrimary()
    }

    // primary := number | s | ( expression )
    private func parsePrimary() throws -> Rational {
        switch peek() {
        case .number(let v):
            advance()
            return try Rational(num: [v], den: [1])
        case .ident(let name):
            advance()
            if name.lowercased() != "s" {
                throw ParseError.message("使える変数は s だけ（今: \(name)）")
            }
            return try Rational(num: [0, 1], den: [1]) // s
        case .lparen:
            advance()
            let r = try parseExpression()
            if peek() != .rparen { throw ParseError.message(") が足りない") }
            advance()
            return r
        case .end:
            throw ParseError.message("有効な式を入力して下さい。")
        default:
            throw ParseError.message("式として読めないトークン: \(peek())")
        }
    }
}

// =====================================================
//  Complex式をそのまま評価する（非整数冪OK）パーサ
//  対応: + - * / ^(実数), 暗黙積, ( ), s, sqrt(x)
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
        case .number, .ident, .lparen:
            return true
        default:
            return false
        }
    }

    func parse() throws -> CExpr {
        let e = try parseExpression()
        if peek() != .end { throw ParseError.message("式の後ろに余計なものがある") }
        return e
    }

    // expression := term { (+|-) term }
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

    // term := power { (*|/|implicitMul) power }
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

    // power := unary { ^ realConst }
    // realConst は 例: 2, -1, 0.5, ( -1/2 ) みたいな「実数定数式」
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

            // sqrt(x) だけサポート（必要なら exp/log も増やせる）
            if lower == "sqrt" {
                guard peek() == .lparen else { throw ParseError.message("sqrt の後ろは ( ) が必要") }
                advance()
                let arg = try parseExpression()
                if peek() != .rparen { throw ParseError.message(") が足りない") }
                advance()
                return .sqrt(arg)
            }

            throw ParseError.message("使える変数は s、関数は sqrt() のみ（今: \(name)）")

        case .lparen:
            advance()
            let e = try parseExpression()
            if peek() != .rparen { throw ParseError.message(") が足りない") }
            advance()
            return e

        case .end:
            throw ParseError.message("有効な式を入力して下さい。")
        default:
            throw ParseError.message("式として読めないトークン: \(peek())")
        }
    }

    // ---- ここから「実数定数式」だけ読む（sは禁止） ----
    private func parseRealConst() throws -> Double {
        if peek() == .lparen {
            advance()
            let v = try parseRealExpr()
            if peek() != .rparen { throw ParseError.message(") が足りない（指数）") }
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
            if peek() != .rparen { throw ParseError.message(") が足りない（指数）") }
            advance()
            return v
        case .ident:
            throw ParseError.message("指数に s は使えない（指数は実数定数だけ）")
        default:
            throw ParseError.message("指数が読めない")
        }
    }
}

// =====================================================
//  UI（式入力 → 係数化 → Bode）
// =====================================================
struct ContentView: View {

    // ✅ 式入力（アプリ再起動しても保持）
    @AppStorage("bode_exprText") private var exprText: String = ""
    //@AppStorage("bode_exprText") private var exprText: String = "G(s)=100/(s*(s+1)^2)"
    @FocusState private var isExprFocused: Bool

    // 周波数レンジ（アプリ再起動しても保持）
    @AppStorage("bode_wStartExp") private var wStartExp: Double = -2
    @AppStorage("bode_wEndExp") private var wEndExp: Double = 2
    @AppStorage("bode_points") private var points: Double = 1000

    // 詳細設定（周波数範囲）を開閉（アプリ再起動しても保持）
    @AppStorage("bode_isAdvancedOpen") private var isAdvancedOpen: Bool = false

    @State private var bodeData: [BodePoint] = []
    @State private var errorMessage: String?

    // 入力変更のたびに再計算（スライダー連打で重くならないよう軽いデバウンス）
    @State private var recalcTask: Task<Void, Never>?

    // デバッグ表示（係数）
    @State private var lastNum: [Double] = []
    @State private var lastDen: [Double] = []

    var body: some View {
        NavigationStack {
            ScrollView {
                VStack(alignment: .leading, spacing: 16) {

                    GroupBox("伝達関数") {
                        VStack(alignment: .leading, spacing: 10) {
                            TextField("G(s)=100/(s*(s+1)^2)", text: $exprText, axis: .vertical)
                                .textFieldStyle(.roundedBorder)
                                .lineLimit(3...6)
                                .focused($isExprFocused)
//                                .overlay(alignment: .trailing) {
//                                    if !exprText.isEmpty {
//                                        Button {
//                                            exprText = ""
//                                            // クリア後にキーボードを閉じたいなら↓を有効化
//                                            // isExprFocused = false
//                                        } label: {
//                                            Image(systemName: "xmark.circle.fill")
//                                                .foregroundStyle(.secondary)
//                                                .padding(.trailing, 10)
//                                        }
//                                        .buttonStyle(.plain)
//                                        .accessibilityLabel("入力をクリア")
//                                    }
//                                }

                            if let errorMessage {
                                Text(errorMessage)
                                    .font(.caption)
                                    .foregroundStyle(.red)
                            }

                            // 係数を見たいなら（不要なら丸ごと消してOK）
                            if !lastNum.isEmpty, !lastDen.isEmpty {
                                HStack(spacing: 12) {
                                    Text("分子: \(lastNum.map { String(format: "%.4g", $0) }.joined(separator: ", "))")
                                    Text("分母: \(lastDen.map { String(format: "%.4g", $0) }.joined(separator: ", "))")
                                }
                                .font(.caption2)
                                .foregroundStyle(.secondary)
                                .lineLimit(1)
                                .truncationMode(.tail)
                            }
                            
//                            Text("points: \(bodeData.count)")
//                                .font(.caption2)
//                                .foregroundStyle(.secondary)
                        }
                    }
                    
                    // リセットボタン
                    Button(role: .destructive) {
                        exprText = ""          // ← 数式をクリア
                        isExprFocused = false  // ← キーボードも閉じる（不要なら消してOK）
                    } label: {
                        Label("リセット", systemImage: "arrow.counterclockwise")
                            .frame(maxWidth: .infinity)
                    }
                    .buttonStyle(.bordered)
                    .disabled(exprText.isEmpty)
                    

                    Button {
                        withAnimation(.easeInOut(duration: 0.2)) {
                            isAdvancedOpen.toggle()
                        }
                    } label: {
                        HStack(spacing: 8) {
                            Image(systemName: isAdvancedOpen ? "chevron.down" : "chevron.right")
                                .font(.caption.weight(.semibold))
                            Text("詳細設定")
                                .font(.body.weight(.semibold))
                            Spacer()
                        }
                        .contentShape(Rectangle())
                    }
                    .buttonStyle(.plain)
                    .padding(.vertical, 4)

                    if isAdvancedOpen {
                        GroupBox("周波数範囲（rad/s）") {
                            VStack(alignment: .leading, spacing: 10) {
                                HStack {
                                    Text("開始 10^")
                                    Slider(value: $wStartExp, in: -6...3, step: 0.5)
                                    Text(String(format: "%.1f", wStartExp))
                                        .monospacedDigit()
                                        .frame(width: 50, alignment: .trailing)
                                }
                                HStack {
                                    Text("終了 10^")
                                    Slider(value: $wEndExp, in: -4...6, step: 0.5)
                                    Text(String(format: "%.1f", wEndExp))
                                        .monospacedDigit()
                                        .frame(width: 50, alignment: .trailing)
                                }
                                HStack {
                                    Text("点数")
                                    Slider(value: $points, in: 200...3000, step: 100)
                                    Text("\(Int(points))")
                                        .monospacedDigit()
                                        .frame(width: 60, alignment: .trailing)
                                }
                            }
                        }
                        .transition(.opacity.combined(with: .move(edge: .top)))
                        
                        Button {
                            resetAdvancedSettings()
                        } label: {
                            Label("設定をリセット", systemImage: "arrow.counterclockwise")
                                .frame(maxWidth: .infinity)
                        }
                        .buttonStyle(.bordered)
                    }

                    BodeChartsView(
                        data: bodeData,
                        wStartExp: min(wStartExp, wEndExp - 0.5),
                        wEndExp: max(wEndExp, wStartExp + 0.5)
                    )
                    
                    Text("※これらの図は片対数グラフであり、横軸は角周波数wを対数目盛(logw)で取っています。")
                        .font(.footnote)
                        .foregroundStyle(.secondary)
                        .padding(.top, 6)
                }
                .padding()
            }
            .navigationTitle("ボード線図")
            .scrollDismissesKeyboard(.interactively)
            .simultaneousGesture(
                TapGesture().onEnded {
                    isExprFocused = false
                }
            )
        }
        .onAppear { scheduleRecalc() }
        .onChange(of: exprText) { _, _ in scheduleRecalc() }
        .onChange(of: wStartExp) { _, _ in scheduleRecalc() }
        .onChange(of: wEndExp) { _, _ in scheduleRecalc() }
        .onChange(of: points) { _, _ in scheduleRecalc() }
    }

    private func scheduleRecalc() {
        recalcTask?.cancel()
        recalcTask = Task {
            // 体感はほぼ即時、でもドラッグ中の連続計算を抑える
            try? await Task.sleep(nanoseconds: 150_000_000)
            if Task.isCancelled { return }
            await MainActor.run {
                recalc()
            }
        }
    }

    private func recalc() {
        errorMessage = nil
        lastNum = []
        lastDen = []

        // start < end を保証
        let start = min(wStartExp, wEndExp - 0.5)
        let end = max(wEndExp, wStartExp + 0.5)

        do {
            // ---- まずは従来の「係数化」ルート（整数冪向け） ----
            let parser = try Parser(exprText)
            let r = try parser.parse()
            let coeffs = r.toHighFirstCoeffs()

            lastNum = coeffs.num
            lastDen = coeffs.den

            let raw = BodeCalc.bode(
                numHighFirst: coeffs.num,
                denHighFirst: coeffs.den,
                wStartExp: start,
                wEndExp: end,
                points: Int(points)
            )

            bodeData = raw.filter { $0.logW.isFinite && $0.magDB.isFinite && $0.phaseDeg.isFinite }
            if bodeData.isEmpty {
                errorMessage = "計算結果が NaN/∞ になってるっぽい（式や周波数範囲を確認してね）"
            }
            return

        } catch {
            // ---- ダメなら「直接評価」ルート（非整数冪OK） ----
            do {
                let p2 = try ComplexExprParser(exprText)
                let expr = try p2.parse()

                lastNum = []
                lastDen = []

                let raw = BodeCalc.bodeEval(
                    wStartExp: start,
                    wEndExp: end,
                    points: Int(points)
                ) { s in
                    expr.eval(s)
                }

                bodeData = raw.filter { $0.logW.isFinite && $0.magDB.isFinite && $0.phaseDeg.isFinite }
                if bodeData.isEmpty {
                    errorMessage = "計算結果が NaN/∞ になっています（式や周波数範囲を確認して下さい）"
                }
            } catch {
                bodeData = []
                errorMessage = (error as? LocalizedError)?.errorDescription ?? "\(error)"
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
//  Charts（ゲイン/位相）
// xは log10(ω) にして擬似log軸
// =====================================================
struct BodeChartsView: View {
    let data: [BodePoint]
    let wStartExp: Double
    let wEndExp: Double
    
    var sorted: [BodePoint] { data.sorted { $0.logW < $1.logW } }
    
    // --- PHASE Y-AXIS HELPERS ---
    private var phaseMinTick: Double {
        let minV = sorted.map { $0.phaseDeg }.min() ?? -180
        return floor(minV / 90.0) * 90.0
    }
    
    private var phaseMaxTick: Double {
        let maxV = sorted.map { $0.phaseDeg }.max() ?? 180
        return ceil(maxV / 90.0) * 90.0
    }
    
    var body: some View {
        VStack(alignment: .leading, spacing: 18) {
            
            GroupBox("ゲイン線図　x:w, y:20log|G(jw)|[dB]") {
                Chart(sorted) { p in
                    LineMark(
                        x: .value("logW", p.logW),
                        y: .value("Mag(dB)", p.magDB)
                    )
                    .interpolationMethod(.catmullRom)
                }
                .frame(height: 240)
                .chartXScale(domain: wStartExp...wEndExp)
                .chartXAxis { bodeXAxis(startExp: wStartExp, endExp: wEndExp) }
                .chartYAxis { AxisMarks(position: .leading) }
            }
            
            GroupBox("位相線図　x:w, y:∠G(jw)[deg]") {
                Chart(sorted) { p in
                    LineMark(
                        x: .value("logW", p.logW),
                        y: .value("Phase(deg)", p.phaseDeg)
                    )
                    .interpolationMethod(.catmullRom)
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
                                Text(String(format: "%.0f", y))
                            }
                        }
                    }
                }
            }
        }
    }
    
    @AxisContentBuilder
    private func bodeXAxis(startExp: Double, endExp: Double) -> some AxisContent {
        let step = 1.0

        // domainの内側に入る整数だけ（ResultBuilder対策でクロージャ化）
        let exps: [Double] = {
            let start = ceil(startExp)
            let end = floor(endExp)

            if start > end {
                // 1 decade未満で整数が入らないときは中央に1個だけ
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
                    Text(format10PowLabel(exp: x))
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
        return "10" + superscriptInt(e)   // 10⁻¹ みたいになる
    }
}


