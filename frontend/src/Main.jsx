import React from 'react';
import { Route, Routes} from 'react-router-dom';


import Home from './pages/Home';
import Tools from './pages/Tools';
import About from './pages/About';

export default function Main() {
    return (
        <Routes>
            <Route exact path='/' element={<Home />}></Route>
            <Route exact path='/tools' element={<Tools />}></Route>
            <Route exact path='/about' element={<About />}></Route>
        </Routes>
    );
}