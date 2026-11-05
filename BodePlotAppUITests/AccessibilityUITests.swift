import XCTest

final class AccessibilityUITests: XCTestCase {
    @MainActor
    func testExpressionAndSettingsAccessibility() throws {
        let app = XCUIApplication()
        app.launchArguments = ["-AppleLanguages", "(en)", "-AppleLocale", "en_US"]
        app.launch()

        // Use the app's controls to start with one empty expression. Run this
        // test on an isolated simulator because it changes persisted data.
        let clear = app.buttons["Clear All"]
        XCTAssertTrue(clear.waitForExistence(timeout: 10))
        if clear.isEnabled {
            for _ in 0..<5 where !clear.isHittable { app.swipeUp() }
            XCTAssertTrue(clear.isHittable)
            clear.tap()
        }

        let add = app.buttons["Add Expression"]
        for _ in 0..<5 where !add.isHittable { app.swipeDown() }
        XCTAssertTrue(add.waitForExistence(timeout: 10))
        add.tap()
        add.tap()

        for index in 1...3 {
            XCTAssertTrue(app.textFields["Transfer Functions \(index)"].waitForExistence(timeout: 10))
            XCTAssertTrue(app.buttons["Delete G\(index)"].exists)
        }

        let settings = app.buttons["Advanced Settings"]
        for _ in 0..<5 where !settings.isHittable { app.swipeUp() }
        XCTAssertTrue(settings.waitForExistence(timeout: 10))
        if settings.value as? String == "Expanded" { settings.tap() }
        XCTAssertEqual(settings.value as? String, "Collapsed")
        settings.tap()
        XCTAssertEqual(settings.value as? String, "Expanded")

        let start = app.sliders["Start 10^"]
        let end = app.sliders["End 10^"]
        let points = app.sliders["Points"]
        XCTAssertTrue(start.waitForExistence(timeout: 10))
        XCTAssertTrue(end.exists)
        XCTAssertTrue(points.exists)

        let reset = app.buttons["Reset to Defaults"]
        for _ in 0..<5 where !reset.isHittable { app.swipeUp() }
        XCTAssertTrue(reset.waitForExistence(timeout: 10))
        reset.tap()
        XCTAssertEqual(start.value as? String, "-2.0")
        XCTAssertEqual(end.value as? String, "2.0")
        XCTAssertEqual(points.value as? String, "1,000")

        if !start.isHittable { app.swipeDown() }
        start.adjust(toNormalizedSliderPosition: 0.5)
        XCTAssertNotEqual(start.value as? String, "-2.0")
    }
}
