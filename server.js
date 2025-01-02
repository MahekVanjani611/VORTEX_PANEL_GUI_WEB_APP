const express = require('express');
const bodyParser = require('body-parser');
const { exec } = require('child_process');
const fs = require('fs');
const path = require('path');

const app = express();
const PORT = 3000;

app.use(bodyParser.json());
app.use(express.static(path.join(__dirname)));

// Serve the main page (input form)
app.get('/', (req, res) => {
    res.sendFile(path.join(__dirname, 'index.html'));
});

// Serve page for streamlines plot
app.get('/streamlines', (req, res) => {
    res.sendFile(path.join(__dirname, 'streamlines.html'));
});

// Endpoint to process inputs and invoke C++ executable
app.post('/calculate', (req, res) => {
    const { ymc, xmc, tm, trailing_edge_type, alfa } = req.body;

    console.log('Received form data:', { ymc, xmc, tm, trailing_edge_type, alfa }); // Log received data

    const command = `.\\aa.exe ${ymc} ${xmc} ${tm} ${trailing_edge_type} ${alfa}`;
    console.log('Executing command:', command); // Log command execution

    exec(command, (error, stdout, stderr) => {
        if (error) {
            console.error(`Error executing the program: ${error.message}`);
            return res.status(500).json({ error: 'Failed to execute program' });
        }

        if (stderr) {
            console.error(`stderr: ${stderr}`);
            return res.status(500).json({ error: 'Error in program execution' });
        }

        console.log('Program executed successfully, reading output files...');

        Promise.all([
            fs.promises.readFile('Cp.csv', 'utf8'),
            fs.promises.readFile('panel_points.csv', 'utf8'),
            fs.promises.readFile('results.csv', 'utf8'),
        ]).then(([cpData, panelData, resultsData]) => {
            console.log('Files read successfully'); // Confirm files are read

            // Parse Cp.csv file
            const cpPoints = cpData.split('\n').map(line => {
                const [x, y] = line.trim().split(/\s+/).map(Number);
                return { x, y: -y }; // Reverse y values for Cp
            }).filter(point => !isNaN(point.x) && !isNaN(point.y));
            console.log('Parsed Cp points:', cpPoints); // Log parsed Cp points

            // Parse panel_points.csv file
            const panelPoints = panelData.split('\n').map(line => {
                const [x, y] = line.trim().split(/\s+/).map(Number);
                return { x, y };
            }).filter(point => !isNaN(point.x) && !isNaN(point.y));
            console.log('Parsed Panel points:', panelPoints); // Log parsed panel points

            // Parse results.csv file
            const resultsLines = resultsData.split('\n').filter(line => line.trim() !== '');
            const results = {};
            resultsLines.forEach(line => {
                const [key, value] = line.split('=');
                results[key.trim()] = parseFloat(value.trim()) || value.trim(); // Store values as float or string
            });
            console.log('Parsed results:', results); // Log parsed results

            // Read and parse streamline files
            const streamlineFiles = [
                'forward_stagline.csv',
                'aft_stagline.csv',
                'panel_points.csv',
                'other_streamlines_0.csv', 'other_streamlines_1.csv', 'other_streamlines_2.csv',
                'other_streamlines_3.csv', 'other_streamlines_4.csv', 'other_streamlines_5.csv',
                'other_streamlines_6.csv', 'other_streamlines_7.csv', 'other_streamlines_8.csv', 'other_streamlines_9.csv'
            ];

            // Read each file in parallel
            const filePromises = streamlineFiles.map(fileName =>
                fs.promises.readFile(fileName, 'utf8').catch(err => null) // Catch errors for missing files
            );

            Promise.all(filePromises).then(fileContents => {
                let streamlinesPoints = [];

                fileContents.forEach((content, index) => {
                    if (content) {
                        // Parse the file with x and y values only
                        const streamlineData = content.split('\n').map(line => {
                            const parts = line.trim().split(/\s+/).map(Number);
                            if (parts.length === 2) { // Handle files with only x and y values
                                const [x, y] = parts;
                                return { x, y };
                            }
                            return null; 
                        }).filter(point => point !== null); // Filter out invalid lines
                        
                        if (streamlineData.length > 0) {
                            streamlinesPoints.push({ file: streamlineFiles[index], data: streamlineData });
                            console.log(`Parsed data for ${streamlineFiles[index]}`, streamlineData); // Log parsed streamline data
                        } else {
                            console.log(`No valid data found in ${streamlineFiles[index]}`);
                        }
                    } else {
                        console.log(`File ${streamlineFiles[index]} not found or is empty.`);
                    }
                });

                const responseData = {
                    cpPoints,
                    panelPoints,
                    results,
                    concatenatedValues: `${ymc}${xmc}${tm}`,
                    writealfa: `${alfa}`,
                    streamlinesPoints // Include parsed streamlines data
                };

                res.json(responseData);
            }).catch(err => {
                console.error(`Error reading streamline files: ${err.message}`);
                res.status(500).json({ error: 'Failed to read streamline files' });
            });

        }).catch(err => {
            console.error(`Error reading files: ${err.message}`);
            res.status(500).json({ error: 'Failed to read output files' });
        });
    });
});

// Start the server
app.listen(PORT, () => {
    console.log(`Server is running at http://localhost:${PORT}`);
});
