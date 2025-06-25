import React from 'react';
import { Route, Routes} from 'react-router-dom';


import Home from './pages/Home';
import Tools from './pages/Tools';
import About from './pages/About';
import ToolPage from "./pages/ToolPage";
import {ScrollRestoration} from "react-router-dom";

export default function Main({toolsList, setToolsList}) {
    return (
        <Routes>
            <Route exact path='/' element={<Home toolsList={toolsList} setToolsList={setToolsList}/>}></Route>
            <Route exact path='/tools' element={<Tools toolsList={toolsList} setToolsList={setToolsList}/>}></Route>
            <Route exact path='/about' element={<About />}></Route>
            <Route exact path='/tool-page/:toolURL' element={<ToolPage toolsList={toolsList} setToolsList={setToolsList}/>}></Route>
        </Routes>
    );
}