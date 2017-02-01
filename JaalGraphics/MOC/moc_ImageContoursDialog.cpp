/****************************************************************************
** Meta object code from reading C++ file 'ImageContoursDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "ImageContoursDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ImageContoursDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JImageContoursDialog_t {
    QByteArrayData data[15];
    char stringdata0[174];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JImageContoursDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JImageContoursDialog_t qt_meta_stringdata_JImageContoursDialog = {
    {
        QT_MOC_LITERAL(0, 0, 20), // "JImageContoursDialog"
        QT_MOC_LITERAL(1, 21, 8), // "readFile"
        QT_MOC_LITERAL(2, 30, 0), // ""
        QT_MOC_LITERAL(3, 31, 10), // "newContour"
        QT_MOC_LITERAL(4, 42, 13), // "setNodeRadius"
        QT_MOC_LITERAL(5, 56, 12), // "setEdgeWidth"
        QT_MOC_LITERAL(6, 69, 12), // "setNodeColor"
        QT_MOC_LITERAL(7, 82, 12), // "setEdgeColor"
        QT_MOC_LITERAL(8, 95, 9), // "deleteAll"
        QT_MOC_LITERAL(9, 105, 10), // "displayIDs"
        QT_MOC_LITERAL(10, 116, 13), // "deleteSegment"
        QT_MOC_LITERAL(11, 130, 12), // "closeContour"
        QT_MOC_LITERAL(12, 143, 6), // "saveAs"
        QT_MOC_LITERAL(13, 150, 11), // "closeDialog"
        QT_MOC_LITERAL(14, 162, 11) // "genEdgeMesh"

    },
    "JImageContoursDialog\0readFile\0\0"
    "newContour\0setNodeRadius\0setEdgeWidth\0"
    "setNodeColor\0setEdgeColor\0deleteAll\0"
    "displayIDs\0deleteSegment\0closeContour\0"
    "saveAs\0closeDialog\0genEdgeMesh"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JImageContoursDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    13,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   79,    2, 0x08 /* Private */,
    3,    0,   80,    2, 0x08 /* Private */,
    4,    0,   81,    2, 0x08 /* Private */,
    5,    0,   82,    2, 0x08 /* Private */,
    6,    0,   83,    2, 0x08 /* Private */,
    7,    0,   84,    2, 0x08 /* Private */,
    8,    0,   85,    2, 0x08 /* Private */,
    9,    0,   86,    2, 0x08 /* Private */,
    10,    0,   87,    2, 0x08 /* Private */,
    11,    0,   88,    2, 0x08 /* Private */,
    12,    0,   89,    2, 0x08 /* Private */,
    13,    0,   90,    2, 0x08 /* Private */,
    14,    0,   91,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JImageContoursDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JImageContoursDialog *_t = static_cast<JImageContoursDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->readFile();
            break;
        case 1:
            _t->newContour();
            break;
        case 2:
            _t->setNodeRadius();
            break;
        case 3:
            _t->setEdgeWidth();
            break;
        case 4:
            _t->setNodeColor();
            break;
        case 5:
            _t->setEdgeColor();
            break;
        case 6:
            _t->deleteAll();
            break;
        case 7:
            _t->displayIDs();
            break;
        case 8:
            _t->deleteSegment();
            break;
        case 9:
            _t->closeContour();
            break;
        case 10:
            _t->saveAs();
            break;
        case 11:
            _t->closeDialog();
            break;
        case 12:
            _t->genEdgeMesh();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JImageContoursDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JImageContoursDialog.data,
        qt_meta_data_JImageContoursDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JImageContoursDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JImageContoursDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JImageContoursDialog.stringdata0))
        return static_cast<void*>(const_cast< JImageContoursDialog*>(this));
    if (!strcmp(_clname, "Ui::ImageContoursDialog"))
        return static_cast< Ui::ImageContoursDialog*>(const_cast< JImageContoursDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JImageContoursDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 13)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 13;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 13)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 13;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
