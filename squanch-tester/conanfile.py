from conans import ConanFile, CMake
import os

class BlitzTestConan(ConanFile):
    name = "squanch-tester"
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"
    requires = "catch2/2.1.1@bincrafters/stable", "blitz/1.0.1@tuncb/pangea", "eigen/3.3.4@conan/stable"

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def imports(self):
        self.copy("*.dll", dst="bin", src="bin")
        self.copy("*.dylib*", dst="bin", src="lib")
        self.copy('*.so*', dst='bin', src='lib')

    def test(self):
        os.chdir("bin")
        self.run(".%sexample" % os.sep)
