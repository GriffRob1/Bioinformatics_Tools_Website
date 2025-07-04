import React from 'react';
import ReactDOM from 'react-dom/client';
import './styles/index.css';
import "@radix-ui/themes/styles.css";
import App from './App';
import reportWebVitals from './testing/reportWebVitals';
import { BrowserRouter } from 'react-router-dom';
import {Theme} from "@radix-ui/themes";


const root = ReactDOM.createRoot(document.getElementById("root"));
root.render(
    <React.StrictMode>
        <BrowserRouter>
            <Theme >
                <App />
            </Theme>
        </BrowserRouter>
    </React.StrictMode>
);

// If you want to start measuring performance in your app, pass a function
// to log results (for example: reportWebVitals(console.log))
// or send to an analytics endpoint. Learn more: https://bit.ly/CRA-vitals
reportWebVitals();
