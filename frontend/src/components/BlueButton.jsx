import React from 'react'
import {useState} from "react";
import {useNavigate} from "react-router-dom";

export default function BlueButton({children, URL, buttonClass, onClick}) {

    const [isMouseOver, setIsMouseover] = useState(false);

    const setBorderWhite = () => {
        setIsMouseover(true);
    }

    const setBorderBlack = () => {
        setIsMouseover(false);
    }

    const buttonStyle = {
        borderColor: isMouseOver ? 'white' : 'black'
    }

    const navigate = useNavigate();
    const onButtonClick = () => {
        if (onClick) {
            onClick();
        }
        navigate(URL);
    }

    return (
        <button className={`blue-button ${buttonClass}`}
                onClick={onButtonClick}
                onMouseOver={setBorderWhite}
                onMouseOut={setBorderBlack}
                style={buttonStyle}
                >
            {children}
        </button>
    )

}