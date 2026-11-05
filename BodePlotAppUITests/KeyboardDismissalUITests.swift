import XCTest

final class KeyboardDismissalUITests: XCTestCase {
    @MainActor
    func testOutsideTapDismissesKeyboardWithoutBlockingControls() {
        let app = XCUIApplication()
        app.launchArguments = ["-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()

        let clear = app.buttons["Clear All"]
        XCTAssertTrue(clear.waitForExistence(timeout: 10))
        if clear.isEnabled { clear.tap() }
        let add = app.buttons["Add Expression"]
        add.tap()

        let first = app.textFields["Transfer Functions 1"]
        let second = app.textFields["Transfer Functions 2"]
        if !second.waitForExistence(timeout: 3) { add.tap() }
        XCTAssertTrue(second.waitForExistence(timeout: 5))
        first.tap()
        XCTAssertTrue(app.keyboards.firstMatch.waitForExistence(timeout: 5))
        second.tap()
        second.typeText("1")
        XCTAssertTrue(app.keyboards.firstMatch.exists)
        XCTAssertEqual(second.value as? String, "1")

        let background = app.coordinate(withNormalizedOffset: .zero).withOffset(
            CGVector(dx: app.frame.midX, dy: first.frame.minY - 25)
        )
        background.tap()
        XCTAssertTrue(app.keyboards.firstMatch.waitForNonExistence(timeout: 5))

        first.tap()
        XCTAssertTrue(app.keyboards.firstMatch.waitForExistence(timeout: 5))
        app.navigationBars["Bode Plot"].tap()
        XCTAssertTrue(app.keyboards.firstMatch.waitForNonExistence(timeout: 5))

        first.tap()
        XCTAssertTrue(app.keyboards.firstMatch.waitForExistence(timeout: 5))
        add.tap()
        XCTAssertTrue(app.textFields["Transfer Functions 3"].waitForExistence(timeout: 5))
        XCTAssertTrue(app.keyboards.firstMatch.waitForNonExistence(timeout: 5))
    }
}
