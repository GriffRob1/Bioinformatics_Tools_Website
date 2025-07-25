import '../styles/NavBar.css';
import React from 'react';
import {NavLink} from 'react-router-dom';
import BlueButton from "./BlueButton";

export default function NavBar() {
//TODO change link to sign in page
    return (
        <nav className={'container nav-bar'}>
            <NavLink to='/'><h1>Bioinformatics Tools</h1></NavLink>
            <ul>
                <li><NavLink to='/'>Home</NavLink></li>
                <li><NavLink to='/tools'>Tools</NavLink></li>
                <li><NavLink to='/about'>About</NavLink></li>
                <li><BlueButton URL={'/'} buttonClass={'sign-in-button'}>Sign In</BlueButton></li>
            </ul>
        </nav>
    )
}