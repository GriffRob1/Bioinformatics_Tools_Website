import React from 'react';
import { NavLink } from 'react-router-dom';

export default function NavBar() {
    return (
        <nav className='container nav-bar'>
            <h1>Bioinformatics Tools</h1>
            <ul>
                <li><NavLink to='/'>Home</NavLink></li>
                <li><NavLink to='/tools'>Tools</NavLink></li>
                <li><NavLink to='/about'>About</NavLink></li>
                <li><button><NavLink to='/'>Sign In</NavLink></button></li>
            </ul>
        </nav>
    )
}