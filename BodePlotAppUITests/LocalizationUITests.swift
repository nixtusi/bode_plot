import XCTest

final class LocalizationUITests: XCTestCase {
    @MainActor
    func testRepresentativeLanguages() throws {
        let cases: [(language: String, locale: String, settings: String, frequency: String, error: String)] = [
            ("en", "en_US", "Advanced Settings", "Frequency Range (rad/s)", "Invalid character: @"),
            ("ja", "ja_JP", "詳細設定", "周波数範囲 (rad/s)", "不正な文字: @"),
            ("de", "de_DE", "Erweiterte Einstellungen", "Frequenzbereich (rad/s)", "Ungültiges Zeichen: @"),
            ("ar", "ar_SA", "الإعدادات المتقدمة", "نطاق التردد (راديان/ثانية)", "حرف غير صالح: @")
        ]

        for item in cases {
            let app = XCUIApplication()
            app.launchArguments = ["-AppleLanguages", "(\(item.language))", "-AppleLocale", item.locale]
            app.launch()

            let settings = app.buttons[item.settings]
            if !settings.isHittable { app.swipeUp() }
            XCTAssertTrue(settings.waitForExistence(timeout: 10), "Settings missing in \(item.language)")
            let frequency = app.staticTexts[item.frequency]
            if !frequency.exists { settings.tap() }
            XCTAssertTrue(frequency.waitForExistence(timeout: 5), "Frequency label missing in \(item.language)")
            if !frequency.isHittable { app.swipeUp() }
            attachScreenshot("\(item.language)-settings")

            app.swipeDown()
            let field = app.textFields.firstMatch
            XCTAssertTrue(field.waitForExistence(timeout: 5), "Expression field missing in \(item.language)")
            field.tap()
            field.typeText("@")
            // Localized interpolation may contain bidirectional isolation marks.
            let errorPrefix = String(item.error.dropLast(3))
            let error = app.staticTexts.matching(NSPredicate(format: "label CONTAINS %@", errorPrefix)).firstMatch
            XCTAssertTrue(error.waitForExistence(timeout: 5), "Parse error missing in \(item.language)")
            XCTAssertTrue(error.label.contains("@"), "Invalid character missing in \(item.language): \(error.label.debugDescription)")
            attachScreenshot("\(item.language)-error")
            app.terminate()
        }
    }

    private func attachScreenshot(_ name: String) {
        let attachment = XCTAttachment(screenshot: XCUIScreen.main.screenshot())
        attachment.name = name
        attachment.lifetime = .keepAlways
        add(attachment)
    }
}
