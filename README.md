
<p align="center">
   <img width="1245" alt="logo_project" src="https://github.com/sgrunber/Project-ppchem-tools-kit/raw/main/docs/source/_static/images/logo_project.png">
   
 </p>


<div align="center">
  <h1 style="font-size: 3em;">🛠️ Tools Kit 🛠 </h1>
</div>


<p align="center">
  <a href="https://github.com/sgrunber/Project-ppchem-tools-kit/actions">
    <img src="https://img.shields.io/badge/build-passing-brightgreen.svg?style=for-the-badge&logo=github&logoColor=white" alt="Build Status">
  </a>
  <a href="https://github.com/sgrunber/Project-ppchem-tools-kit">
    <img src="https://img.shields.io/badge/coverage-95%25-brightgreen.svg?style=for-the-badge&logo=github&logoColor=white" alt="Coverage">
  </a>
  <a href="https://github.com/sgrunber/Project-ppchem-tools-kit/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue.svg?style=for-the-badge&logo=github&logoColor=white" alt="License">
  </a>
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

 <a href="https://github.com/sgrunber/Project-ppchem-tools-kit/tree/main"><img src="https://img.shields.io/badge/GitHub--181717.svg?logo=github&logoColor=white&style=for-the-badge&logoWidth=40&logoColor=white" alt="GitHub"></a>

Next, clone the repository using the following command:

#### a) Using HTTPS (default):

```bash
git clone https://github.com/sgrunber/Project-ppchem-tools-kit.git
```
#### b) Using SSH (recommended for users who have set up SSH):

To clone with SSH, use the following command:

```bash
git clone git@github.com:sgrunber/Project-ppchem-tools-kit.git
```


### 3. ⬇️ Install with pip


After cloning, navigate into the cloned directory and install the package locally using pip:

```bash
cd path/to/Project-ppchem-tools-kit
pip install .
```
This will install the package and its dependencies from the cloned repository, automatically fetching the latest version:



  <a href="https://pypi.org/project/Project-ppchem-tools-kit/"><img src="https://img.shields.io/badge/PyPI-v1.0.0-blue.svg?style=for-the-badge&logo=visual-studio-code&logoColor=white" alt="PyPI"></a>

Alternatively, you can install Tools Kit using pip, which will fetch and install the latest updated version from PyPI:

```bash
pip install Project_ppchem-tools-kit
```

## 🚀 Usage
<a id="usage"></a>

To use Tools Kit and access its features, follow these steps after cloning the repository:

### Using Jupyter Notebook:

<a href="https://jupyter.org/"><img src="https://img.shields.io/badge/Jupyter-Notebook-F37626.svg?style=for-the-badge&logo=visual-studio-code&logoColor=white" alt="Jupyter Notebook"></a>

1. Launch Jupyter Notebook by navigating to the cloned repository directory in your terminal and running the command:
   ```bash
   jupyter notebook
   ```
2. In the Jupyter Notebook interface, navigate to the `notebooks` directory.

3. Open the notebook named `project_report.ipynb`.

4. Follow the instructions inside the notebook to interact with the Tools Kit interface and utilize its features for molecule analysis, graph plotting, error calculation, and more.
This notebook provides an interactive environment for convenient usage of Tools Kit directly within Visual Studio Code.


### Usage with Visual Studio Code:

 <a href="https://code.visualstudio.com/"><img src="https://img.shields.io/badge/Visual%20Studio%20Code-007ACC.svg?style=for-the-badge&logo=visual-studio-code&logoColor=white" alt="VSCode"></a>

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
        - Options: The graph allows customization of various parameters such as changing scales, axis labels, title, displaying maximum values, enabling gridlines, adjusting line types and colors, modifying grid and label axes, and background color.

4. 📊 Graph Maker
   - Input: Excel document with x and y values arranged in two columns
   - Output: Graph plotting all values
      - Options: The graph allows customization of various parameters such as changing scales, axis labels, title, displaying maximum values, enabling gridlines, adjusting line types and colors, modifying grid and label axes, and background color.


5. 🧮 Error Propagation
   - Input: Variables, their values, and their uncertainties.
   - Output: Mean value, its uncertainty, and the result can be copied as LaTeX code for easy inclusion in documents.

6. 🔬 Show Molecule
   - Input: SMILES
   - Output: Molecular structure



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

  <a href="http://www.bibtex.org/"><img src="https://img.shields.io/badge/BibTeX-0.99d-blue.svg?style=for-the-badge&logo=visual-studio-code&logoColor=white" alt="BibTeX"></a>

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

<a href="https://www.chicagomanualofstyle.org/home.html"><img src="https://img.shields.io/badge/Chicago%20Style-17th%20Ed.-5579A1.svg?style=for-the-badge&logo=visual-studio-code&logoColor=white" alt="Chicago Style"></a>

```bash
Méloé Enzinger and Sébastien Grunberg. 2024. Project-ppchem-tools-kit. Computer software. GitHub. https://github.com/username/Project-ppchem-tools-kit.` 
```
## 📜 License

<a href="https://github.com/sgrunber/Project-ppchem-tools-kit/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg?style=for-the-badge&logo=visual-studio-code&logoColor=white" alt="License"></a>

<a id="license"></a>

This project is licensed under the MIT License. See the LICENSE file for details.

## 📸 Screenshots
<a id="screenshots"></a>

**Note**: Graphs have much better quality when saved to your computer.

## Interface
<p align="center">
<img width="700" alt="Capture d’écran 2024-05-18 à 16 44 09" src="https://github.com/sgrunber/Project-ppchem-tools-kit/raw/main/docs/source/_static/images/screen.png">
   </p>
   
## Graph
<p align="center">
<img width="700" alt="graph" src="https://github.com/sgrunber/Project-ppchem-tools-kit/assets/160881864/8d1f8352-3815-4d5d-9545-0a7138f132a5">
    </p>
    
## Linear Regression
<p align="center">
<img width="700" alt="linear_reg" src="https://github.com/sgrunber/Project-ppchem-tools-kit/assets/160881864/eb1a56fb-d801-457b-9a1d-56fc258a6478">
 </p>

## Graph & Linear Regression Settings
<p align="center">
<img width="220" alt="scale" src="https://github.com/sgrunber/Project-ppchem-tools-kit/assets/160881864/8bf14d9a-a573-43fd-bfce-721dccda3415">
<img width="178" alt="graph_option" src="https://github.com/sgrunber/Project-ppchem-tools-kit/assets/160881864/2d08b793-6903-418b-a83d-ae3d24d36171">
<img width="220" alt="Label" src="https://github.com/sgrunber/Project-ppchem-tools-kit/assets/160881864/ad9da31f-5d76-495b-8f81-b3324c2944a6">
 </p>

## Show Molecule
<p align="center">
<img width="398" alt="show_molecule" src="https://github.com/sgrunber/Project-ppchem-tools-kit/assets/160881864/691fda4c-1f1d-4aa3-90ef-88cf759a4f70">
 </p>

