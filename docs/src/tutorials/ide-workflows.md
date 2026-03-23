# IDE Workflows for Julia

> [!TIP]
> Working in a good editor significantly improves the Julia experience. Choose whichever option fits your existing workflow best.

While the Julia REPL is powerful on its own, most users pair it with an IDE that provides syntax highlighting, inline documentation, and — most importantly — tight REPL integration so that code can be sent from a script to a live Julia REPL session with a single keypress combination.

## VS Code/VSCodium

[Visual Studio Code](https://code.visualstudio.com/) is probably the most widely used editor for Julia development.
[VSCodium](https://vscodium.com/) is a fully open-source build of VS Code without Microsoft telemetry, and is compatible with the same Julia extensions.

Install the **[Julia extension](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia)** (`julialang.language-julia`) to unlock:

- Integrated REPL — send lines or selections with `Shift+Enter`
- Inline evaluation results and variable display
- Debugger, profiler, and plot pane
- Workspace variable explorer and data viewer

The extension can be installed from within the IDE via the **Extensions panel** (`Ctrl+Shift+X`) by searching for `Julia`.

:::tip Windows users — shell mode fix
If Julia's shell mode (`;` in the REPL) doesn't work inside VS Code, installing [Git Bash](https://git-scm.com/install/) and adding the following to your VS Code `settings.json` (open via `Ctrl+Shift+P` → *Preferences: Open User Settings (JSON)*) has been reported to fix the issue. This was found via trial and error and may not be the only/best solution:

```json
{
    "terminal.integrated.env.windows": {
        "SHELL": "C:\\Program Files\\Git\\bin\\bash.exe",
        "Path": "${env:Path};C:\\Program Files\\Git\\usr\\bin"
    },
    "terminal.integrated.profiles.windows": {
        "Git Bash": {
            "path": "C:\\Program Files\\Git\\bin\\bash.exe"
        }
    },
    "terminal.integrated.defaultProfile.windows": "Git Bash"
}
```

:::

## Positron

[Positron](https://positron.posit.co/) is an open-source fork of VS Code by Posit (makers of RStudio) with first-class support for both R and Python. It ships a dedicated *Connections*, *Variables*, and *Plots* pane out of the box, making it feel closer to RStudio for researchers coming from an R background.

There is also a dedicated **[positron-julia extension](https://open-vsx.org/extension/ntluong95/positron-julia)** available on Open VSX, installable via the Extensions panel. The `positron-julia` extension is in very early development and may be unstable.

## JetBrains IDEs

[JetBrains](https://www.jetbrains.com/) offer Julia support via the **[Flexible Julia plugin](https://plugins.jetbrains.com/plugin/29356-flexible-julia)**. It should work in any JetBrains IDE (e.g., IntelliJ IDEA, PyCharm) and provides:

- Full language server support (completion, go-to-definition, refactoring)
- Integrated Julia REPL within the IDE
- Package manager UI

> [!NOTE]
> JetBrains IDEs are commercial software, though free community editions and academic licences are available.

## Neovim + Iron.nvim

For more terminal-focused workflows, [Neovim](https://neovim.io/) combined with **[Iron.nvim](https://github.com/Vigemus/iron.nvim)** provides a lightweight but highly capable environment:

- `iron.nvim` opens a Julia REPL split and lets you send lines, visual selections, or entire files with configurable keymaps
- Pair with **[nvim-lspconfig](https://github.com/neovim/nvim-lspconfig)** + the Julia language server (`julials`) for completion and diagnostics
