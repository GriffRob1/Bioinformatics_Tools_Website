import React from 'react';
import {NavLink, useNavigate} from 'react-router-dom';
import BlueButton from "./BlueButton";

export default function NavBar() {

    return (
        <nav className={'container nav-bar'}>
            <NavLink to='/'><h1>Bioinformatics Tools</h1></NavLink>
            <ul>
                <li><NavLink to='/'>Home</NavLink></li>
                <li><NavLink to='/tools'>Tools</NavLink></li>
                <li><NavLink to='/about'>About</NavLink></li>
                <li>
                    <BlueButton content={'Sign In'}
                                URL={'/'} //TODO change link to sign in page
                                buttonClass={'sign-in-button'}/>
                </li>
            </ul>
        </nav>
    )
}