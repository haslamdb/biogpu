// Minimal JSON stub for compilation
// In production, use nlohmann/json or similar library

#ifndef JSON_HPP
#define JSON_HPP

#include <map>
#include <string>
#include <vector>
#include <stdexcept>

namespace nlohmann {

class json {
public:
    json() = default;
    
    static json parse(const std::string& str) {
        // Stub implementation
        throw std::runtime_error("JSON parsing not implemented in stub");
    }
    
    bool contains(const std::string& key) const {
        return false;
    }
    
    json& operator[](const std::string& key) {
        static json dummy;
        return dummy;
    }
    
    const json& operator[](const std::string& key) const {
        static json dummy;
        return dummy;
    }
    
    template<typename T>
    operator T() const {
        return T{};
    }
    
    std::string dump(int indent = -1) const {
        return "{}";
    }
    
    class exception : public std::exception {
    public:
        const char* what() const noexcept override {
            return "JSON exception";
        }
    };
};

} // namespace nlohmann

#endif // JSON_HPP