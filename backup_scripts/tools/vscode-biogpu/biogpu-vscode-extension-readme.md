# BioGPU VS Code Extension Installation Guide

## Quick Installation

### Method 1: Install from VSIX file (Recommended for now)

1. **Create the extension directory structure:**
```bash
cd ~/Documents/Code/biogpu/tools/vscode-biogpu 

```

2. **Save the files:**
- Save `package.json` (from the first artifact) in the root directory
- Save `biogpu.tmLanguage.json` (from the second artifact) in `syntaxes/`
- Save `language-configuration.json` (from the third artifact) in the root
- Save `biogpu-snippets.json` (from the fourth artifact) in `snippets/`

3. **Package the extension:**
```bash
# Install vsce (Visual Studio Code Extension manager)
npm install -g vsce

# Package the extension
vsce package

# This creates: biogpu-language-0.1.0.vsix
```

4. **Install in VS Code:**
- Open VS Code
- Press `Ctrl+Shift+P` (or `Cmd+Shift+P` on Mac)
- Type "Install from VSIX"
- Select the `.vsix` file you created
- Reload VS Code

### Method 2: Local Development Installation

1. **Copy the extension folder to VS Code extensions directory:**

**Windows:**
```bash
copy vscode-biogpu-extension %USERPROFILE%\.vscode\extensions\biogpu-language-0.1.0
```

**macOS/Linux:**
```bash
cp -r vscode-biogpu-extension ~/.vscode/extensions/biogpu-language-0.1.0
```

2. **Reload VS Code**

## Features

### Syntax Highlighting
The extension provides syntax highlighting for:
- Pipeline and stage declarations
- GPU kernel decorators (`@gpu_kernel`, `@parallel`)
- Built-in functions (`parallel_map`, `scan_mutations`, etc.)
- Data types (`fastq_file`, `genome_database`, etc.)
- Keywords and control flow
- Comments and strings

### Code Snippets
Type these prefixes and press Tab:
- `pipeline` - Create a new pipeline
- `stage` - Add a pipeline stage
- `fqresist` - Fluoroquinolone resistance check
- `pmap` - Parallel mapping
- `filter` - Read filtering
- `abundance` - Calculate abundance
- `report` - Clinical report
- `gpu` - GPU kernel decorator

### Auto-completion
- Automatic bracket/quote closing
- Smart indentation
- Code folding for pipelines and stages

## Testing the Extension

Create a test file `test.biogpu`:

```biogpu
# Test BioGPU syntax highlighting
pipeline TestFQResistance {
    input: {
        reads: fastq_file,
        patient_id: string
    }
    
    output: {
        report: json,
        abundance: csv
    }
    
    @gpu_kernel(memory: 16GB)
    stage analyze {
        # This should be highlighted as a comment
        filtered = filter_reads(reads) {
            min_quality: 20,
            min_length: 100
        }
        
        mutations = scan_mutations(filtered, fq_database) {
            positions: {
                gyrA: [83, 87],
                parC: [80, 84]
            }
        }
        
        emit: mutations
    }
}
```

## Customization

### Adding Custom Snippets
Edit `snippets/biogpu-snippets.json` to add your own snippets.

### Modifying Syntax Colors
You can customize colors in your VS Code settings:

```json
{
    "editor.tokenColorCustomizations": {
        "textMateRules": [
            {
                "scope": "keyword.control.biogpu",
                "settings": {
                    "foreground": "#FF6B6B"
                }
            },
            {
                "scope": "entity.name.decorator.biogpu",
                "settings": {
                    "foreground": "#4ECDC4",
                    "fontStyle": "italic"
                }
            }
        ]
    }
}
```

## Troubleshooting

1. **Extension not working?**
   - Make sure file has `.biogpu` extension
   - Try reloading VS Code (`Ctrl+Shift+P` → "Reload Window")

2. **Syntax highlighting looks wrong?**
   - Check that no other extensions are conflicting
   - Verify the file is recognized as "BioGPU" in the status bar

3. **Snippets not appearing?**
   - Make sure you're in a `.biogpu` file
   - Try typing the exact prefix and pressing Tab

## Future Plans

- IntelliSense support for BioGPU functions
- Linting and error checking
- Integration with BioGPU compiler
- Hover documentation
- Go to definition support

## Contributing

To modify the extension:
1. Edit the grammar in `syntaxes/biogpu.tmLanguage.json`
2. Add new snippets in `snippets/biogpu-snippets.json`
3. Test changes by pressing F5 in VS Code (opens new window with extension loaded)

Happy coding with BioGPU! 🧬⚡