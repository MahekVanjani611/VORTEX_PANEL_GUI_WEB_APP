# VORTEX PANEL GUI WEB APP

This project is a web application that integrates with a C++ backend to perform vortex panel method simulations, visualize aerodynamic data, and plot airfoil geometry and streamline distributions. It provides an interactive interface for users to input parameters and visualize results of airfoil simulations.

## Features

- **Airfoil Geometry Plotting**: Visualize the airfoil shape and panel points.
- **Cp Distribution Plot**: Plot the coefficient of pressure (Cp) along the airfoil surface.
- **Streamline Plotting**: Plot streamline data to analyze the flow around the airfoil.
- **C++ Backend Integration**: Communicates with a C++ backend (`a`) for vortex panel simulations and data handling.
- **User Interface**: A clean and responsive web interface to interact with simulation data.

## Prerequisites

To run this project locally, ensure you have the following installed:

- **Node.js** (v14 or higher) - [Install Node.js](https://nodejs.org/)
- **npm** (comes with Node.js)
- **C++ Compiler**: A working C++ compiler (e.g., `g++` on Linux/Mac, MinGW on Windows).
- **Git** - To clone the repository.

## Installation

### 1. Clone the Repository

Clone this repository to your local machine:

```bash
git clone https://github.com/MahekVanjani611/VORTEX_PANEL_GUI_WEB_APP.git
cd VORTEX_PANEL_GUI_WEB_APP
```
# Install Dependencies
Install the required Node.js packages by running:
```bash
npm install
```
# Set Up the C++ Backend
The project relies on a C++ backend (a) for performing vortex panel simulations. Make sure the C++ backend is built and ready to use:

1)If you're on Linux or macOS, use g++ to compile the C++ code:
```bash
g++ -o aa 12nov.cpp
```
2)On Windows, ensure that you have a C++ compiler (e.g., MinGW) set up correctly.
 
# Configuration
Make sure the Node.js server can locate and communicate with the C++ executable (a). Update the path if necessary in the Node.js configuration files like server.js.

# Running the Project Locally
1. Start the C++ Backend
In one terminal window, navigate to the project directory and run the C++ executable:
```bash
.\aa
```
This will start the C++ backend, which the Node.js server will communicate with.

2. Start the Node.js Server
In a separate terminal window, start the Node.js server by running:

```bash
node server.js
```
The server will start locally on http://localhost:3000.

3. Access the Application
Open your browser and go to:
```bash
http://localhost:3000
```
The web interface should now be up, allowing you to interact with the application and visualize vortex panel simulation results.
