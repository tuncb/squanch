from conans import ConanFile

class BlitzConan(ConanFile):
    name = "squanch"
    version = "0.1.0"
    license = "https://www.apache.org/licenses/LICENSE-2.0"
    url = "https://github.com/tuncb/squanch"
    description = ("A 3D NURBS manipulation library in modern C++")
    requires = "blitz/1.0.1@tuncb/pangea", "eigen/3.3.4@conan/stable"
    settings = "os", "compiler", "build_type", "arch"

    def package_id(self):
        self.info.header_only()

    def package(self):
        self.copy("*", dst="squanch", src="./include/squanch")

    def package_info(self):
        self.cpp_info.includedirs = ['.']
