
<p align="center">
   <img width="1245" alt="logo_project" src="https://github.com/sgrunber/Project-ppchem-tools-kit/raw/main/docs/source/_static/images/logo_project.png">

 </p>




 <h1 align="center">🧰 Tools Kit</h1>

 <p align="center">
  <a href="https://github.com/sgrunber/Project-ppchem-tools-kit/actions"><img src="https://img.shields.io/badge/build-passing-brightgreen" alt="Build Status"></a>
  <a href="https://github.com/sgrunber/Project-ppchem-tools-kit"><img src="https://img.shields.io/badge/coverage-95%25-brightgreen" alt="Coverage"></a>
  <a href="https://github.com/sgrunber/Project-ppchem-tools-kit/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue" alt="License"></a>
</p>

## 📖 Description
The **Tools Kit** package allows chemistry students to find basic information on molecules and provides easy access to some graphs and calculations required in experimental laboratories.

**🌙 Note**: For an aesthetically pleasing experience, it is recommended to use the dark mode interface if you are on a Mac or any other platform that supports it.


## 📚 Table of Contents
- [Installation](#️-installation)
- [Usage](#️-usage)
- [Features](#️-features)
- [Contributing](#️-contributing)
- [Cite Us](#️-cite-us)
- [License](#️-license)
- [Screenshots](#️-screenshots)

## 🛠️ Installation
<a id="installation"></a>
### 1. 🍴 Fork the Repository
To start, fork the repository to your own GitHub account.  
To do so, navigate to the repository page and click the *"Fork"* button. The repository will then be copied to your account, and you will be able to access it.

### 2. 📥 Clone the Repository
Next, clone the repository using the following command (replace *username* with your GitHub username):

```bash
git clone https://github.com/*username*/Project-ppchem-tools-kit.git
```

### 3. ⬇️ Install with pip
You can also install Tools Kit using pip:

```bash
pip install Project_ppchem-tools-kit
```

## 🚀 Usage
<a id="usage"></a>

To use Tools Kit and access its features, follow these steps after cloning the repository:

### Using Jupyter Notebook:

1. Launch Jupyter Notebook by navigating to the cloned repository directory in your terminal and running the command:
   ```bash
   jupyter notebook
   ```
2. In the Jupyter Notebook interface, navigate to the `notebooks` directory.

3. Open the notebook named `project_report.ipynb`.

4. Follow the instructions inside the notebook to interact with the Tools Kit interface and utilize its features for molecule analysis, graph plotting, error calculation, and more.
This notebook provides an interactive environment for convenient usage of Tools Kit directly within Visual Studio Code.


### Usage with Visual Studio Code:

1. Open Visual Studio Code and navigate to the cloned directory.

2. Install the Jupyter extension if not already installed. You can do this by searching for "Jupyter" in the Extensions view (Ctrl+Shift+X) and installing the "Python" extension pack.

3. Open the `notebooks` directory in Visual Studio Code.

4. Open the notebook named `project_report.ipynb`.

5. Follow the instructions inside the notebook to interact with the Tools Kit interface and utilize its features for molecule analysis, graph plotting, error calculation, and more.
This notebook provides an interactive environment for convenient usage of Tools Kit with Jupyter Notebook.

## 🧪 Features

1. 🧬 Molecule Name
   - Input: Raw formula of the molecule
   - Output: SMILES of the molecule

2. ⚖️ Molecular Weight
   - Input: SMILES of the molecule
   - Output: Corresponding molar mass in g/mol

3. 📈 Linear Regression
   - Input: Excel document (imported by pressing the Browse button)
   - Output: Linear regression graph with the R<sup>2</sup> value

4. 📊 Graph Maker
   - Input: Excel document with x and y values arranged in two columns
   - Output: Graph plotting all values
      - Options: The graph allows customization of various parameters such as changing scales, axis labels, title, displaying maximum values, enabling gridlines, adjusting line types and colors, modifying grid and label axes, and background color.


5. 📉 Error Propagation
   - Input: Variables, their values, and their uncertainties.
   - Output: Mean value, its uncertainty, and the result can be copied as LaTeX code for easy inclusion in documents.

6. 🧪 Show Molecule
   - Input: SMILES
   - Output: Molecular structure

The graphs in sections 3 and 4 can also include input for the title and X, Y axis labels.

## 🤝 Contributing
<a id="contributing"></a>

Contributions are welcome! To contribute to **Tools Kit**, please follow these steps:

1. **Fork** the project to your GitHub account.
2. Create a new branch for your feature or bug fix: 
   ```bash
   git checkout -b feature/your-feature-name
   ```
3. **Create a new branch for your feature or bug fix:**

    ```bash
    git checkout -b feature/your-feature-name
    ```

4. **Make your changes and ensure they adhere to the project's coding conventions and style guidelines.**

5. **Commit your changes with descriptive messages explaining the purpose of your changes:**

    ```bash
    git commit -am 'Add some feature'
    ```

6. **Push your changes to your branch on your forked repository:**

    ```bash
    git push origin feature/your-feature-name
    ```

7. **Once you've pushed your changes, you can open a Pull Request (PR) in the original repository. Provide a clear description of the changes you've made in the PR.**

8. **After reviewing your PR, if everything looks good, it will be merged into the main repository. Congratulations, you've successfully contributed to *Tools Kit*!**

Thank you for your contribution! Your efforts help improve the project for everyone.

## 📚 Cite Us
<a id="cite-us"></a>

If you use Tools Kit in your research work or academic projects, we would appreciate it if you could cite us. You can use the following BibTeX entry:

```bibtex
@misc{project-ppchem-tools-kit,
  title = {Project-ppchem-tools-kit},
  author = {{Méloé Enzinger & Sébastien Grunberg}},
  year = {2024},
  publisher = {GitHub},
  journal = {GitHub Repository},
  howpublished = {\url{https://github.com/username/Project-ppchem-tools-kit}}
}
```

### Chicago Style
```bash
Enzinger, Méloé, and Sébastien Grunberg. 2024. Project-ppchem-tools-kit. Computer software. GitHub. https://github.com/username/Project-ppchem-tools-kit.` 
```
## 📜 License
<a id="license"></a>

This project is licensed under the MIT License. See the LICENSE file for details.

## 📸 Screenshots
<a id="screenshots"></a>
<img width="1002" alt="Capture d’écran 2024-05-18 à 16 44 09" src="https://github.com/sgrunber/Project-ppchem-tools-kit/raw/main/docs/source/_static/images/screen.png">



