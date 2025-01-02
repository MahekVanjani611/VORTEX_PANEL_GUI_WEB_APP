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

    const command = path.join(__dirname, 'aa.exe') + ` ${ymc} ${xmc} ${tm} ${trailing_edge_type} ${alfa}`;
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

        // Ensure output files exist before reading them
        const filePaths = [
            'Cp.csv', 
            'panel_points.csv', 
            'results.csv', 
            'forward_stagline.csv',
            'aft_stagline.csv',
            'other_streamlines_0.csv',
            'other_streamlines_1.csv',
            'other_streamlines_2.csv',
            'other_streamlines_3.csv',
            'other_streamlines_4.csv',
            'other_streamlines_5.csv',
            'other_streamlines_6.csv',
            'other_streamlines_7.csv',
            'other_streamlines_8.csv',
            'other_streamlines_9.csv'
        ];

        const filePromises = filePaths.map(file => 
            fs.promises.readFile(file, 'utf8').catch(err => null)
        );

        Promise.all(filePromises).then(fileContents => {
            const [cpData, panelData, resultsData, ...streamlineData] = fileContents;

            // Check if the required files exist and handle errors
            if (!cpData || !panelData || !resultsData) {
                return res.status(404).json({ error: 'One or more required output files are missing' });
            }

            // Parse Cp.csv file
            const cpPoints = cpData.split('\n').map(line => {
                const [x, y] = line.trim().split(/\s+/).map(Number);
                return { x, y: -y }; // Reverse y values for Cp
            }).filter(point => !isNaN(point.x) && !isNaN(point.y));
            console.log('Parsed Cp points:', cpPoints);

            // Parse panel_points.csv file
            const panelPoints = panelData.split('\n').map(line => {
                const [x, y] = line.trim().split(/\s+/).map(Number);
                return { x, y };
            }).filter(point => !isNaN(point.x) && !isNaN(point.y));
            console.log('Parsed Panel points:', panelPoints);

            // Parse results.csv file
            const resultsLines = resultsData.split('\n').filter(line => line.trim() !== '');
            const results = {};
            resultsLines.forEach(line => {
                const [key, value] = line.split('=');
                results[key.trim()] = parseFloat(value.trim()) || value.trim();
            });
            console.log('Parsed results:', results);

            // Parse streamlines files
            let streamlinesPoints = [];
            streamlineData.forEach((content, index) => {
                if (content) {
                    const streamlineDataParsed = content.split('\n').map(line => {
                        const parts = line.trim().split(/\s+/).map(Number);
                        if (parts.length === 2) {
                            const [x, y] = parts;
                            return { x, y };
                        }
                        return null; 
                    }).filter(point => point !== null);

                    if (streamlineDataParsed.length > 0) {
                        streamlinesPoints.push({ file: filePaths[index + 3], data: streamlineDataParsed });
                        console.log(`Parsed data for ${filePaths[index + 3]}`, streamlineDataParsed);
                    }
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
    });
});

// Start the server
app.listen(PORT, () => {
    console.log(`Server is running at http://localhost:${PORT}`);
});
